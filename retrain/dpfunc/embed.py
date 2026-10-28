"""DPFunc's residue features: ESM-1b (esm1b_t33_650M_UR50S) layer 31, one 1,280-vector per residue, FP32.

    python -m retrain.dpfunc.embed [--datasets original_dpfunc archipelago development test] [--part i --parts n]

Proteins longer than 1,022 residues are cut into 1,022-residue windows overlapping by at least 256 (the last window ends at
the last residue); each residue's vector is the mean over the windows that cover it. One window per forward pass, FP32,
TF32 off. Writes inputs/dpfunc/esm1b/<seq_id[:2]>/<seq_id>.npy (length x 1,280) for every protein that has none yet.
One GPU per part; needs fair-esm (`pip install fair-esm`). The paper's runs embedded about 18 proteins per second on an
H100.

The paper's DPFunc follows its Methods (ESM-1b); the released code's feature script loads ESM2-650M instead. The window
rule for proteins over 1,022 residues is ours, not the authors'.
"""
import argparse
import time

import numpy as np

from retrain.common import add_argument, data_root, load, log
from retrain.dpfunc.common import DIMENSION, esm_file

DATASETS = ('original_dpfunc', 'archipelago', 'development', 'test')
WINDOW, OVERLAP, LAYER = 1022, 256, 31


def window_starts(length, size=WINDOW, overlap=OVERLAP):
    if length <= size:
        return [0]
    starts = list(range(0, length - size + 1, size - overlap))
    if starts[-1] != length - size:
        starts.append(length - size)
    return starts


def representation(model, tokens, layer):
    """fair-esm's forward pass up to `layer` (later layers and the logits are not computed). Copied from the paper's
    extractor; for ESM-1b with layer < 33 it equals model(tokens, repr_layers=[layer])['representations'][layer]."""
    padding = tokens.eq(model.padding_idx)
    x = model.embed_scale * model.embed_tokens(tokens)
    esm1b = getattr(model, 'model_version', None) == 'ESM-1b'
    token_dropout = getattr(model.args, 'token_dropout', False) if esm1b else model.token_dropout
    if token_dropout:
        x.masked_fill_((tokens == model.mask_idx).unsqueeze(-1), 0.)
        observed = (tokens == model.mask_idx).sum(-1).to(x.dtype) / (~padding).sum(-1)
        x = x * (1 - .15 * .8) / (1 - observed)[:, None, None]
    if esm1b:
        x = x + model.embed_positions(tokens)
        if model.emb_layer_norm_before:
            x = model.emb_layer_norm_before(x)
    x = x * (1 - padding.unsqueeze(-1).type_as(x))
    x = x.transpose(0, 1)
    if not padding.any():
        padding = None
    for block in model.layers[:layer]:
        x, _ = block(x, self_attn_padding_mask=padding, need_head_weights=False)
    if layer == len(model.layers):
        x = model.emb_layer_norm_after(x)
    return x.transpose(0, 1)


def residue_features(model, alphabet, sequence):
    """(length, 1,280) FP32: each residue averaged over the windows that cover it."""
    import torch
    converter = alphabet.get_batch_converter()
    values = np.zeros((len(sequence), DIMENSION), dtype=np.float32)
    counts = np.zeros(len(sequence), dtype=np.int16)
    with torch.inference_mode():
        for start in window_starts(len(sequence)):
            window = sequence[start:start + WINDOW]
            _, _, tokens = converter([('protein', window)])
            result = representation(model, tokens.cuda(), LAYER).float().cpu().numpy()
            values[start:start + len(window)] += result[0, 1:len(window) + 1]
            counts[start:start + len(window)] += 1
    values /= counts[:, None]
    if not np.isfinite(values).all():
        raise FloatingPointError('non-finite ESM-1b features')
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--datasets', nargs='+', default=list(DATASETS))
    parser.add_argument('--part', type=int, default=0)
    parser.add_argument('--parts', type=int, default=1)
    add_argument(parser)
    arguments = parser.parse_args()
    root = data_root(arguments.data)
    import esm
    import torch
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    sequences = {}
    for name in arguments.datasets:
        ds = load(root, name, labels=False)
        sequences.update(zip(ds.seq_ids, ds.sequences))
    mine = sorted(sequences)[arguments.part::arguments.parts]
    todo = [h for h in mine if not esm_file(root, h).exists()]
    total = sum(len(sequences[h]) for h in todo)
    log(f'ESM-1b: {len(todo):,} of {len(mine):,} proteins to embed ({total:,} residues)')
    model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
    model = model.cuda().eval().requires_grad_(False)
    started, done = time.monotonic(), 0
    for i, h in enumerate(todo, 1):
        target = esm_file(root, h)
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = target.with_name(f'{h}.tmp.npy')
        np.save(temporary, residue_features(model, alphabet, sequences[h]))
        temporary.replace(target)
        done += len(sequences[h])
        if i % 1000 == 0 or i == len(todo):
            rate = done / (time.monotonic() - started)
            log(f'{i:,}/{len(todo):,} proteins, {rate:,.0f} residues/s, ETA {(total - done) / rate / 3600:.1f} h')


if __name__ == '__main__':
    main()
