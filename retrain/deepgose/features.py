"""DeepGO-SE's input: the ESM2-3B layer-36 embedding of each protein, averaged over its first 1,022 residues, in FP32.

The computation is deepgo2's deepgo/extract_esm.py extract_esm() (commit 6c77c98), unchanged apart from reading
sequences from the datasets instead of a FASTA file: fair-esm esm2_t36_3B_UR50D, FastaBatchedDataset batches of at most
4,096 tokens, sequences truncated to 1,022 residues, the mean of layer 36 over residues 1..min(L, 1022) (BOS and EOS
excluded), computed on the CPU. TF32 is off. The paper's vectors were checked against this function (relative L2 below
1e-5).

    python -m retrain.deepgose.features development test goa hq archipelago original_deepgose

For every dataset named, writes <root>/inputs/deepgose/features/<dataset>.npy (float32, proteins x 2,560, rows in the
dataset's order) and <dataset>.ids.txt (its seq_ids). Each distinct sequence is embedded once. Needs one GPU (the FP32
model alone takes 11.4 GB) and the fair-esm package (`pip install fair-esm`).
"""
import argparse
import time
from pathlib import Path

import numpy as np

from archipelago.io import read_matrix, write_matrix
from retrain.common import add_argument, data_root, load, log

DIM = 2560
LAYER = 36
TRUNCATION = 1022
TOKS_PER_BATCH = 4096


def folder(root):
    return Path(root) / 'inputs' / 'deepgose' / 'features'


def read(root, dataset):
    """The features of a dataset, rows in its order (fails if a protein has no row)."""
    _, x = read_matrix(folder(root) / dataset.name, dataset.seq_ids, mmap=False)
    x = np.asarray(x)
    if x.shape != (len(dataset.seq_ids), DIM) or x.dtype != np.float32:
        raise ValueError(f'{folder(root) / dataset.name}.npy: shape {x.shape} {x.dtype}, expected (n, {DIM}) float32')
    return x


def embed(sequences):
    """{key: float32 vector of 2,560} for {key: sequence}, as deepgo2 extract_esm()."""
    import torch
    from esm import FastaBatchedDataset, pretrained
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    model, alphabet = pretrained.load_model_and_alphabet('esm2_t36_3B_UR50D')
    model = model.eval().to('cuda')
    labels = list(sequences)
    dataset = FastaBatchedDataset(labels, [sequences[k] for k in labels])
    batches = dataset.get_batch_indices(TOKS_PER_BATCH, extra_toks_per_seq=1)
    loader = torch.utils.data.DataLoader(dataset, collate_fn=alphabet.get_batch_converter(TRUNCATION),
                                         batch_sampler=batches)
    out, started = {}, time.monotonic()
    with torch.no_grad():
        for n, (keys, strings, tokens) in enumerate(loader, 1):
            representations = model(tokens.to('cuda', non_blocking=True), repr_layers=[LAYER],
                                    return_contacts=False)['representations']
            hidden = representations[LAYER].to(device='cpu')
            for i, key in enumerate(keys):
                vector = hidden[i, 1:min(TRUNCATION, len(strings[i])) + 1].mean(0).clone().numpy()
                if not np.isfinite(vector).all():
                    raise ValueError(f'non-finite embedding for {key}')
                out[key] = vector
            if n % 500 == 0 or n == len(batches):
                rate = len(out) / (time.monotonic() - started)
                log(f'{len(out):,}/{len(sequences):,} proteins, {rate:.0f}/s, '
                    f'{(len(sequences) - len(out)) / rate / 60:.0f} min left')
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('datasets', nargs='+', help='dataset names under <root>/datasets')
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    sets = [load(root, name, labels=False) for name in a.datasets]
    sequences = {h: s for ds in sets for h, s in zip(ds.seq_ids, ds.sequences)}
    log(f'{len(sequences):,} distinct proteins in {", ".join(a.datasets)}')
    vectors = embed(sequences)
    for ds in sets:
        write_matrix(folder(root) / ds.name, ds.seq_ids, np.stack([vectors[h] for h in ds.seq_ids]).astype(np.float32))
        log(f'wrote {folder(root) / ds.name}.npy')


if __name__ == '__main__':
    main()
