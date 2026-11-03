"""FP16 ESM2 encoder with fused pointwise kernels and bounded CUDA graphs.

Library GEMMs and SDPA are retained. This path is explicitly selected: it never
changes checkpoint, residue windows or pooling, and never falls back silently.
It targets the model-level rotary implementation in the installed Transformers.
"""
from collections import OrderedDict
import time

import torch
import torch.nn.functional as F
import triton
import triton.language as tl
from triton.language.extra.cuda import libdevice


@triton.jit
def _residual_norm(X, R, W, Bias, Y, Z, D: tl.constexpr, EPS: tl.constexpr,
                   ADD: tl.constexpr, BLOCK: tl.constexpr):
    row = tl.program_id(0)
    c = tl.arange(0, BLOCK)
    x = tl.load(X + row * D + c, c < D, 0).to(tl.float32)
    if ADD:
        r = tl.load(R + row * D + c, c < D, 0).to(tl.float32)
        # ESM rounds its residual addition to FP16 before LayerNorm.
        x = (x + r).to(tl.float16).to(tl.float32)
        tl.store(Y + row * D + c, x, c < D)
    mean = tl.sum(x, 0) / D
    centered = tl.where(c < D, x - mean, 0.)
    variance = tl.sum(centered * centered, 0) / D
    w = tl.load(W + c, c < D, 0).to(tl.float32)
    b = tl.load(Bias + c, c < D, 0).to(tl.float32)
    z = centered * tl.rsqrt(variance + EPS) * w + b
    tl.store(Z + row * D + c, z, c < D)


@triton.jit
def _esm_gelu(X, Y, N: tl.constexpr, BLOCK: tl.constexpr):
    i = tl.program_id(0) * BLOCK + tl.arange(0, BLOCK)
    x = tl.load(X + i, i < N, 0).to(tl.float32)
    # Preserve each FP16 rounding in ESM's custom GELU expression. F.gelu is
    # deliberately not substituted. Disable FMA at launch as well.
    half = (x * 0.5).to(tl.float16).to(tl.float32)
    divided = (x / 1.4142135623730951).to(tl.float16).to(tl.float32)
    erf = libdevice.erf(divided).to(tl.float16).to(tl.float32)
    factor = (1.0 + erf).to(tl.float16).to(tl.float32)
    tl.store(Y + i, half * factor, i < N)


@triton.jit
def _rotary(QKV, Cos, Sin, Q, K, V, S: tl.constexpr, H: tl.constexpr,
            D: tl.constexpr, N: tl.constexpr, BLOCK: tl.constexpr):
    i = tl.program_id(0) * BLOCK + tl.arange(0, BLOCK)
    d = i % D
    head = (i // D) % H
    pos = (i // (D * H)) % S
    batch = i // (D * H * S)
    src = (batch * S + pos) * (3 * H * D) + head * D + d
    other = src + tl.where(d < D // 2, D // 2, -D // 2)
    sign = tl.where(d < D // 2, -1., 1.)
    q = tl.load(QKV + src, i < N, 0).to(tl.float32)
    qr = tl.load(QKV + other, i < N, 0).to(tl.float32)
    q = (q * (D ** -0.5)).to(tl.float16).to(tl.float32)
    qr = (qr * (D ** -0.5)).to(tl.float16).to(tl.float32)
    k = tl.load(QKV + src + H * D, i < N, 0).to(tl.float32)
    kr = tl.load(QKV + other + H * D, i < N, 0).to(tl.float32)
    v = tl.load(QKV + src + 2 * H * D, i < N, 0)
    cos = tl.load(Cos + pos * D + d, i < N, 0).to(tl.float32)
    sin = tl.load(Sin + pos * D + d, i < N, 0).to(tl.float32)
    dst = ((batch * H + head) * S + pos) * D + d
    tl.store(Q + dst, q * cos + sign * qr * sin, i < N)
    tl.store(K + dst, k * cos + sign * kr * sin, i < N)
    tl.store(V + dst, v, i < N)


def residual_norm(x, residual, norm):
    y = torch.empty_like(x) if residual is not None else x
    z = torch.empty_like(x)
    _residual_norm[(x.numel() // x.shape[-1],)](
        x, residual if residual is not None else x, norm.weight, norm.bias, y, z,
        x.shape[-1], norm.eps, residual is not None,
        triton.next_power_of_2(x.shape[-1]), num_warps=8, enable_fp_fusion=False)
    return y, z


def esm_gelu(x):
    y = torch.empty_like(x)
    _esm_gelu[(triton.cdiv(x.numel(), 1024),)](
        x, y, x.numel(), 1024, enable_fp_fusion=False)
    return y


class FusedEsm2:
    """Reuse an already-loaded ESM model; combined QKV weights are exact copies.

    Set release_original_qkv only after baseline measurement: it saves duplicate
    weight memory, but the original model can then no longer run attention.
    """
    def __init__(self, model, release_original_qkv=False):
        if (model.training or model.config.is_decoder or
                model.config.position_embedding_type != 'rotary' or
                not hasattr(model, 'rotary_embeddings')):
            raise ValueError('Expected eval-mode encoder with model-level rotary')
        if next(model.parameters()).dtype != torch.float16:
            raise ValueError('FP16 model required')
        self.model = model
        self.qkv = []
        for layer in model.encoder.layer:
            a = layer.attention.self
            self.qkv.append((torch.cat([a.query.weight, a.key.weight, a.value.weight]),
                             torch.cat([a.query.bias, a.key.bias, a.value.bias])))
            if release_original_qkv:
                a.query = a.key = a.value = None
        self.tables = {}

    def prepare_shape(self, length):
        if length not in self.tables:
            if not 2 <= length <= 1024:
                raise ValueError('ESM residue window exceeds 1022')
            position_ids = torch.arange(length, device='cuda').unsqueeze(0)
            dummy = self.model.embeddings.word_embeddings.weight
            self.tables[length] = self.model.rotary_embeddings(dummy, position_ids)

    def states(self, input_ids, attention_mask, unmasked=False, layers=None):
        """Final states after the closing LayerNorm; with layers, a dict of states per layer number instead. Layer L below the
        final one is the residual stream after block L (Hugging Face hidden_states[L])."""
        b, s = input_ids.shape
        cos, sin = self.tables[s]
        captured = {}
        x = self.model.embeddings(input_ids=input_ids, attention_mask=attention_mask)
        _, normalized = residual_norm(x, None, self.model.encoder.layer[0].attention.LayerNorm)
        for j, (layer, (weight, bias)) in enumerate(zip(self.model.encoder.layer, self.qkv)):
            a = layer.attention.self
            h, d = a.num_attention_heads, a.attention_head_size
            mixed = F.linear(normalized, weight, bias)
            q = torch.empty((b, h, s, d), device=x.device, dtype=x.dtype)
            k, v = torch.empty_like(q), torch.empty_like(q)
            _rotary[(triton.cdiv(q.numel(), 1024),)](
                mixed, cos, sin, q, k, v, s, h, d, q.numel(), 1024,
                enable_fp_fusion=False)
            # Explicit mask remains dynamic under graph replay. Key padding is
            # masked; padded queries may compute states but have zero pool weight.
            mask = None if unmasked else attention_mask[:, None, None, :].bool()
            context = F.scaled_dot_product_attention(q, k, v, attn_mask=mask,
                                                     dropout_p=0.0, scale=1.0)
            context = context.transpose(1, 2).contiguous().view(b, s, h * d)
            projected = layer.attention.output.dense(context)
            x, normalized = residual_norm(projected, x, layer.LayerNorm)
            middle = esm_gelu(layer.intermediate.dense(normalized))
            projected = layer.output.dense(middle)
            norm = (self.model.encoder.layer[j + 1].attention.LayerNorm
                    if j + 1 < len(self.qkv) else self.model.encoder.emb_layer_norm_after)
            x, normalized = residual_norm(projected, x, norm)
            if layers is not None and j + 1 in layers and j + 1 < len(self.qkv):
                captured[j + 1] = x  # residual_norm writes a new tensor, so later blocks do not overwrite it
        if layers is None:
            return normalized
        if len(self.qkv) in layers:
            captured[len(self.qkv)] = normalized
        return captured

    def pooled_layers(self, input_ids, attention_mask, weights, layers):
        """Weighted residue sums (FP32) per requested layer, from one forward pass."""
        unmasked = bool(attention_mask.bool().all().item())
        states = self.states(input_ids, attention_mask, unmasked=unmasked, layers=set(layers))
        if set(states) != set(layers):
            raise ValueError(f'Layers not captured: {sorted(set(layers) - set(states))}')
        return {layer: (states[layer].float() * weights.unsqueeze(-1)).sum(1) for layer in layers}

    def __call__(self, input_ids, attention_mask, weights, unmasked=None):
        if unmasked is None:
            # Resolve once per eager forward. Graph preparation resolves this
            # outside capture; masked/unmasked graphs have distinct cache keys.
            unmasked = bool(attention_mask.bool().all().item())
        states = self.states(input_ids, attention_mask, unmasked=unmasked)
        return (states.float() * weights.unsqueeze(-1)).sum(1)


class GraphedEsm2:
    """At most max_graphs private graph pools; LRU eviction is visible in stats.

    Call prepare once per shape outside timing/capture. All replay inputs and
    pooled output use stable addresses. Returned tensor is overwritten on the
    next replay of that shape; copy/consume it before the next call.
    """
    def __init__(self, engine, max_graphs=7):
        if max_graphs < 1:
            raise ValueError('max_graphs must be positive')
        self.engine, self.max_graphs = engine, max_graphs
        self.graphs = OrderedDict()
        self.stats = []

    def prepare(self, input_ids, attention_mask, weights):
        shape = tuple(input_ids.shape)
        unmasked = bool(attention_mask.bool().all().item())
        key = (shape, unmasked)
        if key in self.graphs:
            return
        if len(self.graphs) == self.max_graphs:
            _, old = self.graphs.popitem(last=False)
            del old
            torch.cuda.empty_cache()
        self.engine.prepare_shape(shape[1])
        buffers = [t.clone() for t in (input_ids, attention_mask, weights)]
        started = time.perf_counter()
        stream = torch.cuda.Stream()
        stream.wait_stream(torch.cuda.current_stream())
        with torch.cuda.stream(stream):
            for _ in range(3):
                self.engine(*buffers, unmasked=unmasked)
        torch.cuda.current_stream().wait_stream(stream)
        torch.cuda.synchronize()
        warmed = time.perf_counter()
        graph = torch.cuda.CUDAGraph()
        with torch.cuda.graph(graph):
            output = self.engine(*buffers, unmasked=unmasked)
        torch.cuda.synchronize()
        self.stats.append({'shape':shape, 'unmasked':unmasked, 'warmup_seconds':warmed-started,
                           'capture_seconds':time.perf_counter()-warmed})
        self.graphs[key] = (buffers, graph, output)

    def __call__(self, input_ids, attention_mask, weights):
        shape = tuple(input_ids.shape)
        unmasked = bool(attention_mask.bool().all().item())
        key = (shape, unmasked)
        if key not in self.graphs:
            raise ValueError(f'Graph shape/mask {key} not prepared; capture explicitly')
        buffers, graph, output = self.graphs.pop(key)
        self.graphs[key] = (buffers, graph, output)
        for dst, src in zip(buffers, (input_ids, attention_mask, weights)):
            if dst.dtype != src.dtype:
                raise ValueError('Graph input dtype changed')
            dst.copy_(src)
        graph.replay()
        return output
