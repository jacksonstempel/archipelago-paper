"""DPFunc: the authors' released training data path and per-epoch development check against ours, same model, loss,
optimizer, inputs and GPU.

    python -m retrain.cost.dpfunc [--data DIR] [--subset 6400] [--dev 100] [--dev-only]

Authors' DPFunc_main.py (CSUBioGroup/DPFunc c5faf4df): every training graph unpickled into RAM as a DGL graph with its
residue features, `[(graph, i, label_row)]` fed to dgl GraphDataLoader (batch 64, shuffle, no drop_last), InterPro counts as a
CSR matrix indexed per batch, AdamW 1e-4, FocalLoss, loss.item() every batch; after each epoch `test_performance_gnn_inter`
on the validation set (pure-Python Fmax that walks the ontology for every protein and term). Ours (retrain/dpfunc/train.py):
graphs built per batch from a memory-mapped bank (Graphs), bank pages dropped after each batch, per-epoch development scores
with the focal loss (evaluate).

The ArchipelaGO bank cannot be held in RAM as the authors' loop does (about 1.5 MB of features per protein, over 350 GB), so
both paths are timed on the same random subset of the ArchipelaGO training proteins and the per-batch times are scaled to
the full set (retrain.cost.table). The model (retrain.dpfunc.model, equal to the released one) is used in both arms.
Sequence: ours epoch 1 (cold pages of the bank), ours epochs 2-3, authors epochs 1-3, ours epoch 4. Then the development
check, ours and the authors', on the first --dev development proteins (the paper: 100, and then all 2,500 with --dev-only).

Inputs: retrain.dpfunc's banks and steps for archipelago and development, and the authors' code at
external/dpfunc/DPFunc/ (model_utils, objective). One H100, about 10 minutes for the epochs.
"""
import argparse
import os
import sys
import time

import numpy as np
import torch

from retrain.common import Ontology, add_argument, data_root
from retrain.dpfunc import train as T


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


class TorchIndexable:
    """Compatibility shim only: the released code indexes a scipy CSR matrix with a torch tensor, which current scipy
    rejects (older scipy accepted it). It converts the index to numpy; nothing else changes."""
    def __init__(self, m):
        self.m = m

    def __getitem__(self, i):
        return self.m[i.numpy()]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--subset', type=int, default=6400, help='training proteins timed')
    parser.add_argument('--dev', type=int, default=100, help='development proteins of the development check')
    parser.add_argument('--dev-only', action='store_true', help='only the development check')
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    sys.path.insert(0, str(root / 'external' / 'dpfunc' / 'DPFunc'))
    from DPFunc.model_utils import test_performance_gnn_inter  # the authors' per-epoch validation
    from DPFunc.objective import AverageMeter
    from dgl.dataloading import GraphDataLoader
    from retrain.dpfunc.model import FocalLoss

    rng = np.random.default_rng(0)
    vocabulary = Ontology(root).vocabulary
    T.seed_everything(0)
    torch.set_num_threads(1)
    train = T.Split(root, 'archipelago', len(vocabulary))
    n = len(train.dataset.seq_ids)
    rows = np.sort(rng.choice(n, size=a.subset, replace=False))
    data = T.Graphs(train.bank, train.rows[rows], train.labels[rows])
    domains = train.domains[rows]
    ours_loader = GraphDataLoader(data, batch_size=T.BATCH, drop_last=False, shuffle=True)
    log('subset', a.subset, 'of', n, 'batches', len(ours_loader))
    model = T.new_model(len(vocabulary), 'cuda')
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-4)
    loss_fn = FocalLoss()

    def epoch(name, loader, release):
        model.train()
        meter = AverageMeter()
        torch.cuda.synchronize()
        t0 = time.time()
        for graph, idx, y in loader:
            graph, y = graph.to('cuda'), y.cuda()
            if release:
                train.bank.release()
            logits = model(T.inter_batch(domains, idx.numpy()), graph, graph.ndata['x'])
            loss = loss_fn(logits, y)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            meter.update(loss.item(), len(y))
        torch.cuda.synchronize()
        s = time.time() - t0
        log(f'EPOCH {name} seconds {s:.1f} ms/batch {1000 * s / len(loader):.1f} loss {meter.avg:.4f}')
        return 1000 * s / len(loader)

    if not a.dev_only:
        ours = [epoch('ours cold', ours_loader, True), epoch('ours warm', ours_loader, True),
                epoch('ours warm', ours_loader, True)]
        t = time.time()
        memory = [data[i] for i in range(len(data))]  # the authors' in-RAM list [(graph, i, label_row)]
        log('built the in-RAM graph list (not part of the epoch) seconds', round(time.time() - t, 1))
        theirs_loader = GraphDataLoader(memory, batch_size=64, drop_last=False, shuffle=True)
        theirs = [epoch('authors', theirs_loader, False) for _ in range(3)]
        ours.append(epoch('ours warm', ours_loader, True))
        log('SUMMARY ms/batch ours cold', round(ours[0], 1), 'ours warm', round(float(np.mean(ours[1:])), 1),
            'authors', round(float(np.mean(theirs)), 1))

    # Per-epoch development check, on the first --dev development proteins.
    dev = T.Split(root, 'development', len(vocabulary))
    drows = np.arange(a.dev)
    dev.dataset.seq_ids, dev.rows, dev.domains, dev.labels = (dev.dataset.seq_ids[:a.dev], dev.rows[drows],
                                                              dev.domains[drows], dev.labels[drows])
    t = time.time()
    T.evaluate(model, dev)
    ours_dev = time.time() - t
    log('DEV ours: evaluate seconds', round(ours_dev, 2), 'for', a.dev, 'proteins')
    dmem = [T.Graphs(dev.bank, dev.rows, dev.labels)[i] for i in range(a.dev)]
    idx_goid = dict(enumerate(vocabulary))
    goid_idx = {g: i for i, g in enumerate(vocabulary)}
    labels = dev.labels.toarray().astype(np.float32)
    os.makedirs('data', exist_ok=True)  # the released function writes into ./data
    t = time.time()
    out = test_performance_gnn_inter(model, GraphDataLoader(dmem, batch_size=64, drop_last=False, shuffle=False),
                                     dev.dataset.seq_ids, TorchIndexable(dev.domains), labels, idx_goid, goid_idx, 'mf',
                                     'cuda')
    authors_dev = time.time() - t
    log('DEV authors: test_performance_gnn_inter seconds', round(authors_dev, 2), 'for', a.dev, 'proteins; Fmax',
        round(float(out[0]), 4))
    log('DEV per protein: ours', round(ours_dev / a.dev, 4), 's, authors', round(authors_dev / a.dev, 4), 's')


if __name__ == '__main__':
    main()
