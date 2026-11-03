"""Embed the proteins of one or more datasets with one of the paper's encoders: the extraction that made the paper's
embeddings, unchanged apart from comments and the parts that only read and wrote the paper's files.

    python scripts/embed.py --set archipelago development test --encoder esm2_15b_l48 [--data DIR]

Encoders:
  esm2_15b_l48  ESM2-15B (facebook/esm2_t48_15B_UR50D at revision 5fbca396), layer 48 (after the closing LayerNorm): the
                FP16 encoder with fused pointwise kernels (scripts/esm2_fused_inference.py; matrix products and attention are
                PyTorch's), windows of at most 1,022 residues overlapping by at least 256, each residue averaged over the
                windows covering it, mean over residues in FP32; BOS and EOS excluded
  esm2_3b_l36   ESM2-3B (facebook/esm2_t36_3B_UR50D at revision 476b6399), layer 36, extracted the same way (the recipe
                ablation's input)
  prott5_l21    ProtT5-XL-U50 encoder (Rostlab/prot_t5_xl_half_uniref50-enc at revision 94a6abc0, FP16), the output of
                block 21: residues separated by spaces, U/Z/O/B mapped to X, the first 1,024 tokens, masked mean over every
                attended token including the closing special token
  prott5_l24    the same, its last layer (after the closing LayerNorm; the variant's input)

Work as in the paper: the distinct proteins of the given sets, longest first (ties by seq_id), in units of 1,000 (ESM2) or
2,000 (ProtT5) proteins; in a unit, ESM2 windows are sorted by length and run 8 per forward pass, ProtT5 sequences in
length-sorted batches of at most 8,192 padded tokens. Padding depends on which proteins share a batch, so the last FP16
digits depend on the sets embedded together: the paper embedded HQ, ArchipelaGO's picks, development and test together
(--set archipelago development test) and, separately, each further set (hq_random's picks, the control and dose sets'
added proteins, the published methods' original data). Embedded together as the paper did, whole work units came out
bit for bit equal to the paper's vectors; embedded in other batches they differ slightly (ESM2-15B: cosine 0.9999996 for
32 proteins re-embedded in the paper's own check; ProtT5: the development set embedded on its own, cosine 0.99997, 445
of 2,500 vectors identical). The paper's environments: PyTorch 2.10.0 (CUDA 12.6) and Triton 3.6.0;
Transformers 5.9.0 for ESM2 (the fused encoder needs its model-level rotary embedding) and 4.57.6 with sentencepiece
0.1.99 for ProtT5. Before any work the fused encoder is checked against the Transformers forward on three proteins.

Reads datasets/<set>/sequences.fasta under the data root (--data or $ARCHIPELAGO_DATA). Writes, per set,
inputs/embeddings/<encoder>/<set>.npy (FP16, proteins x dimension, rows in dataset order) and <set>.ids.txt (seq_ids),
the .ids.txt last. One GPU (ESM2-15B: an 80 GB GPU).
"""
import argparse
import math
import re
import sys
import time
from pathlib import Path

import numpy as np

from retrain.common import load
from retrain.paths import data_root

ESM2_15B = dict(model='facebook/esm2_t48_15B_UR50D', revision='5fbca39631164edc1d402a5aa369f982f72ee282')
ESM2_3B = dict(model='facebook/esm2_t36_3B_UR50D', revision='476b639933c8baad5ad09a60ac1a87f987b656fc')
PROTT5 = dict(model_id='Rostlab/prot_t5_xl_half_uniref50-enc', revision='94a6abc029ae13029317b140b7424e012bf8dfbf',
              final_layer=24, dimension=1024)
CONFIG = {
    'esm2_15b_l48': dict(backend='esm2', layer=48, dimension=5120, task_proteins=1000, batch_windows=8, **ESM2_15B),
    'esm2_3b_l36': dict(backend='esm2', layer=36, dimension=2560, task_proteins=1000, batch_windows=8, **ESM2_3B),
    'prott5_l21': dict(backend='t5', layer=21, dimension=1024, task_proteins=2000, token_budget=8192),
    'prott5_l24': dict(backend='t5', layer=24, dimension=1024, task_proteins=2000, token_budget=8192),
}
SUPPORTED = frozenset('ACDEFGHIKLMNPQRSTVWYXBUZO')


def task_order(ids, lengths):
    """Longest proteins first, ties by identifier, so work units have similar lengths."""
    return sorted(range(len(ids)), key=lambda i: (-lengths[i], ids[i]))


def compare(vectors, reference):
    """Largest absolute and mean absolute difference, and smallest and median cosine, over protein pairs."""
    vectors, reference = np.asarray(vectors, dtype=np.float32), np.asarray(reference, dtype=np.float32)
    difference = np.abs(vectors - reference)
    cosine = (vectors * reference).sum(1) / (np.linalg.norm(vectors, axis=1) * np.linalg.norm(reference, axis=1))
    return dict(max_abs=float(difference.max()), mean_abs=float(difference.mean(axis=1).max()), cosine=float(cosine.min()),
                median_cosine=float(np.median(cosine)))


def agree(metrics, cosine=0.9999):
    """Backend self-checks compare two computations of the same states: they must point the same way."""
    if metrics['cosine'] < cosine:
        raise ValueError(f'Backend check failed: {metrics}')
    return metrics


# ----------------------------------------------------------------------------- ESM2

def windows(length):
    if length < 1:
        raise ValueError('Empty sequence')
    if length <= 1022:
        return [(0, length)]
    count = math.ceil((length - 1022) / 766) + 1
    return [(i * (length - 1022) // (count - 1), 1022) for i in range(count)]


def window_weights(length):
    spans_ = windows(length)
    coverage = np.zeros(length, dtype=np.int32)
    for start, size in spans_:
        coverage[start:start + size] += 1
    assert np.all(coverage > 0)
    return [(start, (1.0 / coverage[start:start + size]).astype(np.float32)) for start, size in spans_]


class Esm2:
    def __init__(self, c):
        import torch
        from transformers import AutoModel, AutoTokenizer
        from esm2_fused_inference import FusedEsm2
        torch.set_num_threads(1)
        self.tokenizer = AutoTokenizer.from_pretrained(c['model'], revision=c['revision'])
        self.encoder = AutoModel.from_pretrained(c['model'], revision=c['revision'], dtype=torch.float16,
                                                 attn_implementation='sdpa', add_pooling_layer=False,
                                                 low_cpu_mem_usage=True).eval().cuda()
        with torch.inference_mode():
            self.engine = FusedEsm2(self.encoder)  # the original attention weights stay until release(), for check()
        self.layers, self.dimension, self.batch = [c['layer']], c['dimension'], c['batch_windows']
        if not 1 <= c['layer'] <= self.encoder.config.num_hidden_layers:
            raise ValueError(f'Layer outside 1-{self.encoder.config.num_hidden_layers}')

    def embed(self, sequences):
        import torch
        records, pending = [], []
        for row, sequence in enumerate(sequences):
            weights = window_weights(len(sequence))
            pending.append(len(weights))
            records.extend((row, sequence[start:start + len(w)], w) for start, w in weights)
        records.sort(key=lambda record: len(record[1]))
        totals = {layer: np.zeros((len(sequences), self.dimension), dtype=np.float32) for layer in self.layers}
        with torch.inference_mode():
            for offset in range(0, len(records), self.batch):
                batch = records[offset:offset + self.batch]
                encoded = self.tokenizer([record[1] for record in batch], return_tensors='pt', padding=True, add_special_tokens=True)
                weights = np.zeros(tuple(encoded['input_ids'].shape), dtype=np.float32)
                for j, (_, window, w) in enumerate(batch):
                    if int(encoded['attention_mask'][j].sum()) != len(window) + 2:
                        raise ValueError('Tokenizer changed the residue count')
                    weights[j, 1:len(window) + 1] = w
                tokens = encoded['input_ids'].cuda()
                self.engine.prepare_shape(tokens.shape[1])
                values = self.engine.pooled_layers(tokens, encoded['attention_mask'].cuda(), torch.from_numpy(weights).cuda(), self.layers)
                for layer, value in values.items():
                    value = value.cpu().numpy()
                    for j, (row, _, _) in enumerate(batch):
                        totals[layer][row] += value[j]
                for row, _, _ in batch:
                    pending[row] -= 1
        if any(pending):
            raise ValueError('Windows left unembedded')
        lengths = np.array([len(x) for x in sequences], dtype=np.float32)[:, None]
        return {layer: total / lengths for layer, total in totals.items()}

    def check(self, sequences):
        """The layer against the Hugging Face forward on up to three single-window proteins (before release)."""
        import torch
        single = [x for x in sequences if len(x) <= 1022][:3]
        fused = self.embed(single)
        found = {layer: [] for layer in self.layers}
        with torch.inference_mode():
            for sequence in single:
                encoded = self.tokenizer([sequence], return_tensors='pt', add_special_tokens=True)
                out = self.encoder(input_ids=encoded['input_ids'].cuda(), attention_mask=encoded['attention_mask'].cuda(),
                                   output_hidden_states=True)
                if len(out.hidden_states) != self.encoder.config.num_hidden_layers + 1:
                    raise ValueError('Unexpected number of hidden states')
                for layer in self.layers:
                    found[layer].append(out.hidden_states[layer][0, 1:len(sequence) + 1].float().mean(0).cpu().numpy())
        return {str(layer): agree(compare(fused[layer], np.stack(found[layer]))) for layer in self.layers}

    def release(self):
        """Drop the original attention weights (the fused copies remain); check() cannot run afterwards."""
        import torch
        for layer in self.encoder.encoder.layer:
            attention = layer.attention.self
            attention.query = attention.key = attention.value = None
        torch.cuda.empty_cache()


# ----------------------------------------------------------------------------- ProtT5

def make_batches(records, token_budget):
    records = sorted(records, key=lambda item: (min(len(item[1]) + 2, 1024), item[0]))
    batches, current = [], []
    for record in records:
        proposed = current + [record]
        padded = len(proposed) * min(len(record[1]) + 2, 1024)
        if current and padded > token_budget:
            batches.append(current)
            current = [record]
        else:
            current = proposed
    if current:
        batches.append(current)
    return batches


def pool_attention(hidden, mask):
    import torch
    values = hidden.float()
    weights = mask.to(values.device, dtype=torch.float32).unsqueeze(-1)
    return (values * weights).sum(1) / weights.sum(1).clamp_min(1.0)


class T5LayerCapture:
    def __init__(self, model, layers, capture_final):
        encoder = model.encoder
        if len(encoder.block) < max(layers):
            raise ValueError('Requested layer exceeds loaded encoder depth.')
        self.values = {}
        self.handles = [encoder.block[layer - 1].register_forward_hook(self._hook(layer)) for layer in layers]
        if capture_final:
            self.handles.append(encoder.final_layer_norm.register_forward_hook(self._hook('final')))

    def _hook(self, key):
        def save(_module, _inputs, output):
            self.values[key] = output[0] if isinstance(output, tuple) else output
        return save

    def begin(self):
        self.values = {}

    def finish(self, required):
        missing = set(required).difference(self.values)
        if missing:
            raise RuntimeError(f'Forward pass did not capture: {sorted(map(str, missing))}')
        result, self.values = self.values, {}
        return result


class T5:
    def __init__(self, c):
        import torch
        from transformers import T5EncoderModel, T5Tokenizer
        torch.set_num_threads(1)
        torch.backends.cuda.matmul.allow_tf32 = False
        self.budget = c['token_budget']
        self.tokenizer = T5Tokenizer.from_pretrained(PROTT5['model_id'], revision=PROTT5['revision'], do_lower_case=False)
        self.model = T5EncoderModel.from_pretrained(PROTT5['model_id'], revision=PROTT5['revision'], torch_dtype=torch.float16,
                                                    low_cpu_mem_usage=True).half().eval().cuda()
        if len(self.model.encoder.block) != PROTT5['final_layer'] or c['dimension'] != PROTT5['dimension']:
            raise ValueError('Unexpected ProtT5 depth or width')
        self.layer, self.dimension = c['layer'], c['dimension']
        if not 1 <= self.layer <= PROTT5['final_layer']:
            raise ValueError(f"Layers outside 1-{PROTT5['final_layer']}")
        self.name = 'final' if self.layer == PROTT5['final_layer'] else self.layer
        inner = [self.layer] if self.layer < PROTT5['final_layer'] else []
        self.capture = T5LayerCapture(self.model, inner or [1], capture_final=self.name == 'final')

    def tokenize(self, texts):
        values = [' '.join(re.sub(r'[UZOB]', 'X', text.upper())) for text in texts]
        encoded = self.tokenizer(values, add_special_tokens=True, padding=True, truncation=True, max_length=1024, return_tensors='pt')
        return {key: value.pin_memory().cuda(non_blocking=True) for key, value in encoded.items()}

    def embed(self, sequences):
        import torch
        values = np.zeros((len(sequences), self.dimension), dtype=np.float32)
        for batch in make_batches(list(enumerate(sequences)), self.budget):
            rows, texts = zip(*batch)
            encoded = self.tokenize(list(texts))
            self.capture.begin()
            with torch.inference_mode():
                self.model(**encoded)
            captured = self.capture.finish([self.name])
            values[list(rows)] = pool_attention(captured[self.name], encoded['attention_mask']).cpu().numpy()
        return {self.layer: values}

    def check(self, sequences):
        return {}  # the hooks capture the block outputs that Transformers returns as hidden states

    def release(self):
        pass


BACKENDS = {'esm2': Esm2, 't5': T5}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--set', nargs='+', required=True, help='dataset folder names under datasets/, embedded together')
    parser.add_argument('--encoder', choices=sorted(CONFIG), required=True)
    parser.add_argument('--data', help='data root (default: $ARCHIPELAGO_DATA)')
    args = parser.parse_args()
    root = data_root(args.data)
    c = CONFIG[args.encoder]
    sys.path.insert(0, str(Path(__file__).resolve().parent))  # esm2_fused_inference

    datasets = [load(root, name, labels=False) for name in args.set]
    proteins = {}
    for ds in datasets:
        for h, s in zip(ds.seq_ids, ds.sequences):
            proteins.setdefault(h, s.upper())
    ids = list(proteins)
    sequences = [proteins[h] for h in ids]
    bad = [h for h, s in proteins.items() if not s or set(s) - SUPPORTED]
    if bad:
        raise ValueError(f'{len(bad)} sequences are empty or hold symbols outside {"".join(sorted(SUPPORTED))}, e.g. {bad[0]}')
    order = task_order(ids, [len(s) for s in sequences])
    backend = BACKENDS[c['backend']](c)
    print('CHECK', backend.check(sequences), flush=True)
    backend.release()
    out = root / 'inputs/embeddings' / args.encoder
    out.mkdir(parents=True, exist_ok=True)
    stored = {}  # per set: its matrix (a memory map written as units finish) and each protein's rows in it
    for ds in datasets:
        Path(f'{out / ds.name}.ids.txt').unlink(missing_ok=True)  # written last: its presence marks a complete matrix
        rows = {}
        for r, h in enumerate(ds.seq_ids):
            rows.setdefault(h, []).append(r)
        stored[ds.name] = (np.lib.format.open_memmap(f'{out / ds.name}.npy', mode='w+', dtype=np.float16,
                                                     shape=(len(ds.seq_ids), c['dimension'])), rows)
    started, done = time.monotonic(), 0
    for start in range(0, len(order), c['task_proteins']):
        unit = order[start:start + c['task_proteins']]
        values = backend.embed([sequences[i] for i in unit])[c['layer']].astype(np.float16)
        if not np.isfinite(values).all():
            raise ValueError('Embeddings overflow FP16')
        for matrix, rows in stored.values():
            for j, i in enumerate(unit):
                for r in rows.get(ids[i], ()):
                    matrix[r] = values[j]
        done += len(unit)
        rate = done / max(time.monotonic() - started, 1e-6)
        print(f"{' '.join(args.set)} {args.encoder}: {done:,}/{len(ids):,} proteins, {rate:.1f}/s, "
              f"eta {(len(ids) - done) / rate / 60:.1f} min", flush=True)
    for ds in datasets:
        stored[ds.name][0].flush()
        Path(f'{out / ds.name}.ids.txt').write_text('\n'.join(ds.seq_ids) + '\n', encoding='utf-8')
        print(f'wrote {out / ds.name}.npy: {len(ds.seq_ids):,} proteins', flush=True)


if __name__ == '__main__':
    main()
