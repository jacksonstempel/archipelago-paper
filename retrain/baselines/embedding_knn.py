"""Embedding-kNN baseline: the k training proteins nearest by cosine similarity of ESM2-15B layer-48 mean embeddings; a
term's score is sum_i s_i y_i / sum_i s_i over those neighbours (s_i the similarity, y_i the neighbour's label). k is
chosen from 1, 3, 5, 10, 20, 50 by development micro-AUPRC (ties to the smaller k; the paper's runs chose 20).

    python -m retrain.baselines.embedding_knn --training-sets goa hq archipelago [--device cuda] [--data DIR]

Embeddings: inputs/embeddings/esm2_15b_l48/{<training set>,development,test} (scripts/embed.py). Minutes per training set
on one GPU.
"""
import argparse

import numpy as np
import torch

from archipelago.io import read_matrix
from retrain.baselines import Sets
from retrain.common import add_argument, data_root, log

ENCODER = 'esm2_15b_l48'
K = (1, 3, 5, 10, 20, 50)


def unit(root, name, ids, device):
    _, x = read_matrix(root / 'inputs' / 'embeddings' / ENCODER / name, ids)
    x = torch.as_tensor(np.asarray(x, dtype=np.float32), device=device)
    if not bool(torch.isfinite(x).all()):
        raise ValueError(f'{name}: non-finite embeddings')
    return torch.nn.functional.normalize(x, dim=1)


def neighbours(queries, train, k):
    """Top-k cosine similarities and training rows of every query."""
    sims, rows = [], []
    for start in range(0, len(queries), 512):
        s, r = torch.topk(queries[start:start + 512] @ train.T, k, dim=1)
        sims.append(s)
        rows.append(r)
    return torch.cat(sims), torch.cat(rows)


def knn_scores(sims, rows, labels, k):
    s, r = sims[:, :k], rows[:, :k]
    if not bool((s > 0).all()):
        raise ValueError('a neighbour with non-positive cosine similarity: the weighting is undefined')
    out = torch.empty((len(s), labels.shape[1]), device=s.device)
    for start in range(0, len(s), 256):
        w = s[start:start + 256]
        out[start:start + 256] = (w[:, :, None] * labels[r[start:start + 256]].float()).sum(1) / w.sum(1, keepdim=True)
    return out.cpu().numpy()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--training-sets', nargs='+', required=True)
    parser.add_argument('--device', default='cuda' if torch.cuda.is_available() else 'cpu')
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    for name in a.training_sets:
        sets = Sets(root, name)
        train = unit(root, name, sets.train.seq_ids, a.device)
        labels = torch.as_tensor(sets.train.labels.toarray(), dtype=torch.uint8, device=a.device)
        dev_sims, dev_rows = neighbours(unit(root, 'development', sets.dev.seq_ids, a.device), train, max(K))
        sweep = {k: sets.dev_auprc(knn_scores(dev_sims, dev_rows, labels, k)) for k in K}
        for k, value in sweep.items():
            log(f'embedding_knn {name} k={k}: development micro-AUPRC {value:.4f}')
        chosen = max(K, key=lambda k: (sweep[k], -k))
        test_sims, test_rows = neighbours(unit(root, 'test', sets.test.seq_ids, a.device), train, chosen)
        sets.write('embedding_knn', knn_scores(dev_sims, dev_rows, labels, chosen),
                   knn_scores(test_sims, test_rows, labels, chosen), k=chosen, sweep=sweep)
        log(f'embedding_knn {name}: k = {chosen}')


if __name__ == '__main__':
    main()
