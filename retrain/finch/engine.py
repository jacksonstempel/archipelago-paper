"""The engine that trained Finch's readouts for the paper: one hidden layer on a frozen protein embedding, AdamW with fused
adaptive sharpness-aware minimization (ASAM), each update a CUDA graph with the batch gathered inside it. Needs a CUDA GPU
and Triton.

This is the paper's training code, unchanged apart from comments. Linear layers run in BF16 from shadows of the FP32 master
weights, which the update writes; LayerNorm, the loss, the gradient norms and AdamW are FP32 (AdamW's scalar arithmetic in
FP64, as fused torch AdamW). One update: forward, loss and logit gradient (a Triton kernel), backward; with ASAM, the
gradient norm ||g * s|| (s = |w| + eta), the perturbation rho * s^2 * g / ||g * s|| of the BF16 shadows (LayerNorm
parameters in place), a second forward and backward with the same dropout mask; then gradient-norm clipping at 1 and AdamW
on the FP32 master weights.

Recipe keys (train): `loss` 'bce' (binary cross-entropy on the logits, mean over all protein-term pairs; Finch and every
model of the recipe ablation) or 'asl' (the asymmetric loss, gamma- 4, gamma+ 1, probability margin 0.05; the engine's
default, which the development-set comparison in Methods set aside for binary cross-entropy); `rho` None trains with AdamW alone (one pass,
then the same clipped update); `width` None is the linear readout (LayerNorm and one linear layer to the terms; no hidden
layer, so `dropout` must be None).
"""
import math
import time

import numpy as np
import torch
import torch.nn.functional as F
import triton
import triton.language as tl
from triton.language.extra.cuda import libdevice

MARGIN = 0.05  # ASL probability margin (clip); gamma- 4 and gamma+ 1 are written into _loss_grad


@triton.jit
def _loss_grad(Z, Y, G, L, N: tl.constexpr, VALID_LABELS: tl.constexpr, WIDTH: tl.constexpr, BLOCK: tl.constexpr):
    """ASL (gamma- 4, gamma+ 1, margin 0.05) loss and logit gradient, mean over the N x VALID_LABELS pairs; the output is
    padded to WIDTH columns and padding never contributes."""
    i = tl.program_id(0)*BLOCK + tl.arange(0, BLOCK)
    col = i % WIDTH
    valid = (i < N*WIDTH) & (col < VALID_LABELS)
    j = (i//WIDTH)*VALID_LABELS + col
    z = tl.load(Z+i, i < N*WIDTH, 0).to(tl.float32)
    y = tl.load(Y+j, valid, 0).to(tl.float32)
    p = tl.sigmoid(z)
    neg = tl.sigmoid(-z)
    q = tl.minimum(neg+0.05, 1.)
    logp = tl.minimum(z, 0.) - tl.log(1.+tl.exp(-tl.abs(z)))
    logq = tl.log(tl.maximum(q, 1.e-8))
    a = 1.-p
    b = 1.-q
    b2 = b*b
    b4 = b2*b2
    lp = -a*logp
    ln = -b4*logq
    gp = p*(1.-p)*logp - a*neg
    dq = tl.where(neg+0.05 <= 1., -neg*(1.-neg), 0.)
    gn = 4.*b2*b*dq*logq - b4*tl.where(q >= 1.e-8, dq/tl.maximum(q, 1.e-8), 0.)
    loss = tl.where(valid, tl.where(y > 0., lp, ln), 0.)/(N*VALID_LABELS)
    grad = tl.where(valid, tl.where(y > 0., gp, gn), 0.)/(N*VALID_LABELS)
    tl.store(G+i, grad, i < N*WIDTH)
    tl.store(L+tl.program_id(0), tl.sum(loss, 0))


@triton.jit
def _bce_grad(Z, Y, G, L, N: tl.constexpr, VALID_LABELS: tl.constexpr, WIDTH: tl.constexpr, BLOCK: tl.constexpr):
    """Binary cross-entropy on logits and its logit gradient (sigmoid(z) - y), mean over the N x VALID_LABELS pairs, padded
    like _loss_grad."""
    i = tl.program_id(0)*BLOCK + tl.arange(0, BLOCK)
    col = i % WIDTH
    valid = (i < N*WIDTH) & (col < VALID_LABELS)
    j = (i//WIDTH)*VALID_LABELS + col
    z = tl.load(Z+i, i < N*WIDTH, 0).to(tl.float32)
    y = tl.load(Y+j, valid, 0).to(tl.float32)
    loss = tl.where(valid, tl.maximum(z, 0.) - z*y + tl.log(1.+tl.exp(-tl.abs(z))), 0.)/(N*VALID_LABELS)
    grad = tl.where(valid, tl.sigmoid(z) - y, 0.)/(N*VALID_LABELS)
    tl.store(G+i, grad, i < N*WIDTH)
    tl.store(L+tl.program_id(0), tl.sum(loss, 0))


@triton.jit
def _norm(G, PART, N: tl.constexpr, OFFSET: tl.constexpr, BLOCK: tl.constexpr):
    i = tl.program_id(0)*BLOCK + tl.arange(0, BLOCK)
    g = tl.load(G+i, i < N, 0).to(tl.float32)
    tl.store(PART+OFFSET+tl.program_id(0), tl.sum(g*g, 0))


@triton.jit
def _reduce_norm(PART, NORM, STEP, N: tl.constexpr, BLOCK: tl.constexpr,
                 LOSS=None, VALID=None, CHECK: tl.constexpr = False, ADVANCE: tl.constexpr = True):
    i = tl.arange(0, BLOCK)
    v = tl.load(PART+i, i < N, 0)
    norm = tl.sqrt(tl.sum(v, 0))
    tl.store(NORM, norm)
    if ADVANCE:
        tl.store(STEP, tl.load(STEP)+1.)
    if CHECK:
        loss = tl.load(LOSS)
        tl.store(VALID, (tl.abs(norm) < float('inf')) & (tl.abs(loss) < float('inf')))


@triton.jit
def _asam_norm(G, P, PART, N: tl.constexpr, OFFSET: tl.constexpr, ETA: tl.constexpr, BLOCK: tl.constexpr):
    """Partial squared norm of (abs(parameter) + eta) * gradient."""
    i = tl.program_id(0)*BLOCK + tl.arange(0, BLOCK)
    ok = i < N
    g = tl.load(G+i, ok, 0).to(tl.float32)
    p = tl.load(P+i, ok, 0).to(tl.float32)
    weighted = g*(tl.abs(p)+ETA)
    tl.store(PART+OFFSET+tl.program_id(0), tl.sum(weighted*weighted, 0))


@triton.jit
def _asam_perturb(P, G, ORIGINAL, SHADOW, ADAPTIVE_NORM, N: tl.constexpr,
                  HAS_SHADOW: tl.constexpr, ETA: tl.constexpr, RHO: tl.constexpr, BLOCK: tl.constexpr):
    """Elementwise ASAM step. Linear layers (HAS_SHADOW): only their BF16 forward shadow is perturbed (the second pass reads
    only the shadow) and the FP32 master keeps its value for the update. LayerNorm parameters are read in FP32: saved to
    ORIGINAL and perturbed in place."""
    i = tl.program_id(0)*BLOCK + tl.arange(0, BLOCK)
    ok = i < N
    p = tl.load(P+i, ok, 0).to(tl.float32)
    g = tl.load(G+i, ok, 0).to(tl.float32)
    scale = tl.abs(p)+ETA
    denominator = tl.load(ADAPTIVE_NORM).to(tl.float32)+1.e-12
    if HAS_SHADOW:
        tl.store(SHADOW+i, (p + RHO*g*scale*scale/denominator).to(tl.bfloat16), ok)
    else:
        tl.store(ORIGINAL+i, p, ok)
        tl.store(P+i, p + RHO*g*scale*scale/denominator, ok)


@triton.jit
def _gather_supervision(Y, I, BY, COLS: tl.constexpr, BLOCK: tl.constexpr):
    row = tl.program_id(0)
    col = tl.program_id(1)*BLOCK+tl.arange(0, BLOCK)
    source = tl.load(I+row)
    value = tl.load(Y+source*COLS+col, col < COLS, 0)
    tl.store(BY+row*COLS+col, value, col < COLS)


@triton.jit
def _update(P, G, M, V, SHADOW, NORM, STEP, LR, N: tl.constexpr,
            HAS_SHADOW: tl.constexpr, WD: tl.constexpr, EPS: tl.constexpr, BLOCK: tl.constexpr):
    """AdamW on the FP32 master after global-norm clipping at 1; also writes the BF16 shadow of linear layers."""
    i = tl.program_id(0)*BLOCK + tl.arange(0, BLOCK)
    ok = i < N
    p = tl.load(P+i, ok, 0)
    g = tl.load(G+i, ok, 0).to(tl.float32)
    norm = tl.load(NORM)
    g = g * tl.minimum(1., 1./(norm+1.e-6))
    # Match fused Torch AdamW's double scalar arithmetic, then FP32 storage.
    beta1 = tl.full((), .9, tl.float64)
    beta2 = tl.full((), .999, tl.float64)
    gd = g.to(tl.float64)
    m = (beta1*tl.load(M+i, ok, 0).to(tl.float64) + (1.-beta1)*gd).to(tl.float32)
    v = (beta2*tl.load(V+i, ok, 0).to(tl.float64) + (1.-beta2)*gd*gd).to(tl.float32)
    step = tl.load(STEP)
    lr = tl.load(LR).to(tl.float64)
    b1 = (1.-libdevice.pow(beta1, step.to(tl.float64))).to(tl.float32)
    b2 = tl.sqrt(1.-libdevice.pow(beta2, step.to(tl.float64))).to(tl.float32)
    p = (p.to(tl.float64)-lr*WD*p.to(tl.float64)).to(tl.float32)
    size = (lr/b1.to(tl.float64)).to(tl.float32)
    denom = (tl.div_rn(tl.sqrt(v), b2).to(tl.float64)+EPS).to(tl.float32)
    p = p-tl.div_rn(size*m, denom)
    tl.store(P+i, p, ok)
    tl.store(M+i, m, ok)
    tl.store(V+i, v, ok)
    if HAS_SHADOW:
        tl.store(SHADOW+i, p, ok)


@triton.jit
def _update_original(P, G, M, V, SHADOW, ORIGINAL, NORM, STEP, LR, N: tl.constexpr,
                     HAS_SHADOW: tl.constexpr, WD: tl.constexpr, EPS: tl.constexpr, BLOCK: tl.constexpr):
    """_update starting from ORIGINAL (the unperturbed LayerNorm parameters saved by _asam_perturb)."""
    i = tl.program_id(0)*BLOCK+tl.arange(0, BLOCK)
    ok = i < N
    p = tl.load(ORIGINAL+i, ok, 0).to(tl.float32)
    g = tl.load(G+i, ok, 0).to(tl.float32)
    g = g*tl.minimum(1., 1./(tl.load(NORM)+1.e-6))
    b1 = tl.full((), .9, tl.float64)
    b2 = tl.full((), .999, tl.float64)
    gd = g.to(tl.float64)
    m = (b1*tl.load(M+i, ok, 0).to(tl.float64)+(1.-b1)*gd).to(tl.float32)
    v = (b2*tl.load(V+i, ok, 0).to(tl.float64)+(1.-b2)*gd*gd).to(tl.float32)
    step = tl.load(STEP)
    lr = tl.load(LR).to(tl.float64)
    c1 = (1.-libdevice.pow(b1, step.to(tl.float64))).to(tl.float32)
    c2 = tl.sqrt(1.-libdevice.pow(b2, step.to(tl.float64))).to(tl.float32)
    p = (p.to(tl.float64)-lr*WD*p.to(tl.float64)).to(tl.float32)
    size = (lr/c1.to(tl.float64)).to(tl.float32)
    p = p-tl.div_rn(size*m, (tl.div_rn(tl.sqrt(v), c2).to(tl.float64)+EPS).to(tl.float32))
    tl.store(P+i, p, ok)
    tl.store(M+i, m, ok)
    tl.store(V+i, v, ok)
    if HAS_SHADOW:
        tl.store(SHADOW+i, p, ok)


def loss_gradient(logits, y, loss='asl'):
    n, width = logits.shape
    valid = y.shape[1]
    assert logits.is_contiguous() and y.is_contiguous() and y.shape == (n, valid) and valid <= width
    grad = torch.empty_like(logits)
    partial = torch.empty((triton.cdiv(n*width, 1024),), device=logits.device, dtype=torch.float32)
    {'asl': _loss_grad, 'bce': _bce_grad}[loss][(len(partial),)](logits, y, grad, partial, n, valid, width, 1024)
    return partial.sum(), grad


class Head(torch.nn.Module):
    """LayerNorm, linear, GELU, dropout, linear (archipelago.readout.Readout); hidden_dim None: LayerNorm, linear."""
    def __init__(self, input_dim, hidden_dim, output_dim, dropout):
        super().__init__()
        if hidden_dim is None:
            if dropout is not None:
                raise ValueError('the linear readout has no hidden layer to drop out')
            self.network = torch.nn.Sequential(torch.nn.LayerNorm(input_dim), torch.nn.Linear(input_dim, output_dim))
            return
        self.network = torch.nn.Sequential(
            torch.nn.LayerNorm(input_dim), torch.nn.Linear(input_dim, hidden_dim), torch.nn.GELU(),
            torch.nn.Dropout(dropout), torch.nn.Linear(hidden_dim, output_dim))

    def forward(self, inputs):
        return self.network(inputs)


@torch.no_grad()
def predict(head, x):
    """Sigmoid scores of FP32 inputs in chunks of 256, head in eval mode."""
    previous = head.training
    head.eval()
    result = torch.cat([head(b).sigmoid() for b in x.split(256)])
    head.train(previous)
    return result


class ReadoutStep:
    """One training update (forward, loss, ASAM perturbation, second pass, clipped AdamW) as a CUDA graph. The FP32 head
    holds the master weights; linear layers run in BF16 from shadows written by every update. rho None skips ASAM (one
    pass), loss 'bce' replaces the ASL, and a head without a hidden layer (Head with hidden_dim None) is accepted."""
    def __init__(self, head, lr, rho, eta, weight_decay=0., eps=1e-8, loss='asl'):
        if rho is not None and not (rho > 0 and eta > 0):
            raise ValueError('ASAM rho and eta must be positive')
        if loss not in ('asl', 'bce'):
            raise ValueError(f'unknown loss {loss}')
        self.head, self.wd, self.eps, self.loss = head, float(weight_decay), float(eps), loss
        self.asam = rho is not None
        self.rho, self.eta = (float(rho), float(eta)) if self.asam else (None, None)
        self.lr = torch.as_tensor(lr, device='cuda', dtype=torch.float32)  # filled with the scheduled value every update
        self.params = list(head.parameters())
        if any(p.dtype != torch.float32 for p in self.params):
            raise ValueError('FP32 master parameters required')
        self.linear = len(head.network) == 2  # LayerNorm, linear
        dim = head.network[0].normalized_shape[0]
        self.hidden = None if self.linear else head.network[1].out_features
        self.labels = int(self.params[-1].shape[0])
        self.padded = 128*math.ceil(self.labels/128)
        expected = [(self.labels, dim), (self.labels,)] if self.linear else [(self.hidden, dim), (self.hidden,),
                                                                             (self.labels, self.hidden), (self.labels,)]
        if [tuple(p.shape) for p in self.params][2:] != expected:
            raise ValueError('Unexpected head architecture')
        self.shadows = []
        for p in self.params[2:]:
            shape = (self.padded, *p.shape[1:]) if p.shape[0] == self.labels else p.shape
            shadow = torch.zeros(shape, device=p.device, dtype=torch.bfloat16)
            shadow.flatten()[:p.numel()].copy_(p.detach().flatten())
            shadow.requires_grad_(True)
            self.shadows.append(shadow)
        self.m = [torch.zeros_like(p) for p in self.params]
        self.v = [torch.zeros_like(p) for p in self.params]
        self.blocks = [triton.cdiv(p.numel(), 2048) for p in self.params]
        self.partial = torch.empty(sum(self.blocks), device='cuda')
        self.norm = torch.zeros((), device='cuda')
        self.adaptive_norm = torch.zeros((), device='cuda')
        self.step = torch.zeros((), device='cuda')
        self.valid = torch.empty((), device='cuda', dtype=torch.bool)
        self.originals = [torch.empty_like(p) for p in self.params]
        self.graph = None

    def forward(self, x, dropout_mask):
        x = self.head.network[0](x).to(torch.bfloat16)
        if self.linear:
            return F.linear(x, self.shadows[0], self.shadows[1])
        x = F.gelu(F.linear(x, self.shadows[0], self.shadows[1]))
        return F.linear(x*dropout_mask, self.shadows[2], self.shadows[3])

    def _backward(self, x, y, dropout_mask):
        for p in [*self.params, *self.shadows]:
            p.grad = None
        logits = self.forward(x, dropout_mask)
        loss, gradient = loss_gradient(logits, y, self.loss)
        logits.backward(gradient)
        return loss

    def _gradients(self):
        # Linear gradients stay BF16 (padded rows never enter the norm or the update).
        return [self.params[0].grad, self.params[1].grad, *[s.grad for s in self.shadows]]

    def _reduce(self, kernel_args, target, loss, advance):
        offset = 0
        for p, g, n in zip(self.params, self._gradients(), self.blocks):
            kernel_args(p, g, n, offset)
            offset += n
        _reduce_norm[(1,)](self.partial, target, self.step, sum(self.blocks), triton.next_power_of_2(sum(self.blocks)),
                           loss, self.valid, True, advance)
        torch._assert_async(self.valid, 'Nonfinite loss or gradient norm')

    def _step(self, x, y):
        # Same dropout realization in both passes isolates the weight perturbation.
        dm = None if self.linear else F.dropout(torch.ones((x.shape[0], self.hidden), device=x.device, dtype=torch.bfloat16),
                                                p=self.head.network[3].p, training=True)
        loss = self._backward(x, y, dm)
        if self.asam:
            with torch.no_grad():
                # ASAM: ||g*s|| (s = |w| + eta) sets the radius; displacement rho * s^2 * g / ||g*s||.
                self._reduce(lambda p, g, n, offset: _asam_norm[(n,)](g, p, self.partial, p.numel(), offset, self.eta, 2048),
                             self.adaptive_norm, loss, advance=False)
                for j, (p, g, original) in enumerate(zip(self.params, self._gradients(), self.originals)):
                    shadow = self.shadows[j-2] if j >= 2 else p
                    _asam_perturb[(triton.cdiv(p.numel(), 1024),)](p, g, original, shadow, self.adaptive_norm, p.numel(),
                                                                   j >= 2, self.eta, self.rho, 1024)
            loss = self._backward(x, y, dm)
        self._reduce(lambda p, g, n, offset: _norm[(n,)](g, self.partial, p.numel(), offset, 2048), self.norm, loss,
                     advance=True)
        for j, (p, g, m, v) in enumerate(zip(self.params, self._gradients(), self.m, self.v)):
            grid = (triton.cdiv(p.numel(), 1024),)
            if j < 2 and not self.asam:  # LayerNorm without ASAM: never perturbed
                _update[grid](p, g, m, v, p, self.norm, self.step, self.lr, p.numel(), False, self.wd, self.eps, 1024)
            elif j < 2:  # LayerNorm: perturbed in place; update from the saved original
                _update_original[grid](p, g, m, v, p, self.originals[j], self.norm, self.step, self.lr, p.numel(), False,
                                       self.wd, self.eps, 1024)
            else:
                _update[grid](p, g, m, v, self.shadows[j-2], self.norm, self.step, self.lr, p.numel(), True,
                              self.wd, self.eps, 1024)
        return loss.detach(), self.norm

    def gathered(self, indices, embeddings, targets):
        """One update on rows `indices` of the GPU-resident embeddings (any float dtype, widened to FP32) and dense uint8
        targets, gathered inside the graph."""
        if self.graph is None:
            buffer = indices.clone()
            n, k = buffer.numel(), targets.shape[1]

            def update():
                y = torch.empty((n, k), device=targets.device, dtype=targets.dtype)
                _gather_supervision[(n, triton.cdiv(k, 2048))](targets, buffer, y, k, 2048)
                return self._step(embeddings[buffer].float(), y)
            self.graph = self._capture(update) + (buffer, embeddings.data_ptr(), targets.data_ptr())
        graph, output, buffer, x_ptr, y_ptr = self.graph
        if (x_ptr, y_ptr) != (embeddings.data_ptr(), targets.data_ptr()):
            raise ValueError('The captured graph belongs to other inputs')
        buffer.copy_(indices)
        graph.replay()
        return output

    def _capture(self, update):
        """Record update() as a CUDA graph; the warm-up's changes to weights, moments, step and RNG are undone."""
        values = [*self.params, *self.shadows, *self.m, *self.v, self.step]
        saved = [p.detach().clone() for p in values]
        rng = torch.cuda.get_rng_state()
        stream = torch.cuda.Stream()
        stream.wait_stream(torch.cuda.current_stream())
        with torch.cuda.stream(stream):
            for _ in range(3):
                update()
        torch.cuda.current_stream().wait_stream(stream)
        graph = torch.cuda.CUDAGraph()
        with torch.cuda.graph(graph):
            output = update()
        with torch.no_grad():
            for p, v in zip(values, saved):
                p.copy_(v)
        torch.cuda.set_rng_state(rng)
        return graph, output


class MicroAUPRC:
    """Exact micro-AUPRC (sklearn's non-interpolated AP over all pairs, ties as one threshold group) against fixed labels."""
    def __init__(self, labels):
        labels = np.asarray(labels)
        if labels.ndim != 2 or not np.isin(labels, (0, 1)).all():
            raise ValueError('binary label matrix required')
        self.shape = labels.shape
        self.targets = torch.as_tensor(labels.astype(bool), device='cuda').reshape(-1)
        self.positives = int(self.targets.sum())

    def __call__(self, scores):
        row = torch.as_tensor(scores, dtype=torch.float32, device='cuda').reshape(-1)
        if row.numel() != self.targets.numel() or not bool(torch.isfinite(row).all()):
            raise ValueError('scores must be finite and match the labels')
        positive, counts = torch.unique(row[self.targets], sorted=True, return_counts=True)
        bins = torch.searchsorted(positive, row, right=True)  # bin i+1: scores >= positive[i] and < positive[i+1]
        total = torch.bincount(bins, minlength=len(positive)+1)[1:].flip(0).cumsum(0).flip(0)
        tp = counts.flip(0).cumsum(0).flip(0)
        return float((counts.double()*tp.double()/total.double()).sum()/self.positives)

    def fmax(self, scores, thresholds=np.round(np.arange(1, 101)/100, 2)):
        """Micro Fmax over the evaluation thresholds 0.01-1.00 (predicted positive: score >= threshold)."""
        row = torch.as_tensor(scores, dtype=torch.float32, device='cuda').reshape(-1)
        best = 0.
        for t in thresholds:
            predicted = row >= float(t)
            tp = int((predicted & self.targets).sum())
            best = max(best, 2*tp/(int(predicted.sum())+self.positives))
        return best


class UniformSampler:
    """Training rows drawn uniformly with replacement: the draws of Generator.choice(n, size, p=uniform)."""
    def __init__(self, n):
        self.cdf = np.cumsum(np.full(n, 1/n, dtype=np.float64))
        self.cdf /= self.cdf[-1]

    def draw(self, rng, size):
        return np.searchsorted(self.cdf, rng.random(size), side='right')


def train(x, y, dev_x, dev_y, recipe, seed, log=None):
    """Fit one readout.

    x: training embeddings on the GPU (n, d), any float dtype; y: dense uint8 training labels on the GPU (n, terms); dev_x:
    FP32 development embeddings on the GPU; dev_y: development labels (NumPy, binary). recipe: width, dropout, lr, rho, eta,
    batch, warmup, check_interval, patience, weight_decay, eps, loss (width and dropout None for the linear readout, rho
    None for AdamW alone). Checkpoint: the best raw development micro-AUPRC, checked every
    check_interval updates; training stops after `patience` checks without improvement. Returns (best state dict on the CPU,
    best development scores on the GPU, summary)."""
    r = recipe
    torch.backends.cuda.matmul.allow_tf32 = False
    metric = MicroAUPRC(dev_y)
    sampler = UniformSampler(x.shape[0])
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    rng = np.random.default_rng(seed)
    head = Head(x.shape[1], r['width'], y.shape[1], r['dropout']).cuda()
    engine = ReadoutStep(head, r['lr'], r['rho'], r['eta'], r['weight_decay'], r['eps'], r['loss'])
    best, best_step, stale, step, history = -1., 0, 0, 0, []
    started = time.monotonic()
    torch.cuda.synchronize()
    block = time.monotonic()
    while stale < r['patience']:
        if step % r['check_interval'] == 0:
            # One draw per check interval: the same NumPy draws as one per update.
            drawn = torch.as_tensor(sampler.draw(rng, r['check_interval']*r['batch']).reshape(r['check_interval'], r['batch']),
                                    device='cuda')
        step += 1
        engine.lr.fill_(r['lr']*min(step/r['warmup'], 1))  # linear warm-up, then constant
        loss, norm = engine.gathered(drawn[(step-1) % r['check_interval']], x, y)
        if step % r['check_interval']:
            continue
        torch.cuda.synchronize()
        seconds = time.monotonic()-block
        scores = predict(head, dev_x)
        ap = metric(scores)
        event = dict(update=step, dev_raw_micro_auprc=ap, loss=float(loss), gradient_norm=float(norm),
                     updates_per_second=r['check_interval']/seconds, seconds=time.monotonic()-started)
        if not np.isfinite([ap, event['loss'], event['gradient_norm']]).all():
            raise FloatingPointError(event)
        history.append(event)
        if ap > best:
            # Kept on the GPU, written once after the fit.
            best, best_step, stale, best_scores = ap, step, 0, scores
            best_state = {k: v.detach().clone() for k, v in head.state_dict().items()}
        else:
            stale += 1
        if log:
            log(event, best, best_step)
        block = time.monotonic()
    summary = dict(updates=step, best_update=best_step, dev_raw_micro_auprc=best, seconds=time.monotonic()-started,
                   median_updates_per_second=float(np.median([h['updates_per_second'] for h in history[1:]] or [0])),
                   history=history)
    del engine
    return {k: v.cpu() for k, v in best_state.items()}, best_scores, summary


def load_head(state, dropout=0.):
    """A head rebuilt from a saved state dict (eval mode; dropout plays no part in prediction)."""
    if 'network.4.weight' in state:
        width, dim = state['network.1.weight'].shape
        head = Head(dim, width, state['network.4.weight'].shape[0], dropout).cuda()
    else:  # linear readout
        terms, dim = state['network.1.weight'].shape
        head = Head(dim, None, terms, None).cuda()
    head.load_state_dict(state)
    return head.eval()
