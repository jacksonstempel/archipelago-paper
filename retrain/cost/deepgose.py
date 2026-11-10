"""DeepGO-SE's official training step timed on ArchipelaGO (237,797 proteins) with the axioms of the paper's runs and, in the
same job, with the GO-Plus axioms of an earlier GO release (March 2026; their measured native epoch times on the other
training sets are what the table scales from).

    python -m retrain.cost.deepgose [--data DIR] [--march AXIOMS]

The step is the official one (deepgo2 train.py, commit 6c77c98): DeepGOModel and load_normal_forms (retrain/deepgose/
model.py, unchanged copies), the official FastTensorDataLoader (shuffle = gather of the whole dense feature and label tensors
once per epoch), batch 37, torch Adam lr 5e-4, FP32 (TF32 off): zero_grad, net(x), BCE, net.el_loss(all four normal forms),
backward, step. Two modes on the same model:
  off       losses stay on the GPU; no per-batch syncs
  released  the official train.py bookkeeping: loss.detach().item() and el_loss.detach().item() every batch
The axiom sets are timed in the order paper, march, march, paper, to cancel drift (paper only without --march). Time is
taken over chunks of CHUNK batches with a synchronize at the chunk edges only. Also timed, once: the shuffle, the
development pass and the official per-epoch compute_roc.

Inputs: datasets/{archipelago,development}, their ESM2-3B features (retrain.deepgose.features), the normalized axioms of the
paper's runs (inputs/deepgose/axioms/go-plus.norm, retrain/deepgose/axioms/build.sh) and the authors' code at
external/deepgose/deepgo2/ (its FastTensorDataLoader). --march: the GO-Plus axioms of GO release 2026-03-25, normalized the
same way (not part of the release). One H100, about 10 minutes; needs scikit-learn.
"""
import argparse
import sys
import time

import numpy as np
import torch

from retrain.common import Ontology, add_argument, data_root, load
from retrain.deepgose import features
from retrain.deepgose.model import DeepGOModel, load_normal_forms

WARM, OFF, REL, CHUNK = 10, 300, 150, 50


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def compute_roc(labels, preds):  # deepgo/metrics.py compute_roc, verbatim
    from sklearn.metrics import auc, roc_curve
    fpr, tpr, _ = roc_curve(labels.flatten(), preds.flatten())
    return auc(fpr, tpr)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--march', help='normalized GO-Plus axioms of GO release 2026-03-25 (the earlier native timings)')
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    sys.path.insert(0, str(root / 'external' / 'deepgose' / 'deepgo2'))
    from deepgo.torch_utils import FastTensorDataLoader  # the official loader

    dev = 'cuda:0'
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    axioms = {'paper': str(root / 'inputs' / 'deepgose' / 'axioms' / 'go-plus.norm')}
    if a.march:
        axioms['march'] = a.march
    terms = {t: i for i, t in enumerate(Ontology(root).vocabulary)}
    ds = load(root, 'archipelago', sequences=False)
    n = len(ds.seq_ids)
    log('proteins', n, 'terms', len(terms), 'gpu', torch.cuda.get_device_name(0))
    t = time.time()
    x = torch.from_numpy(np.ascontiguousarray(features.read(root, ds), dtype=np.float32))
    y = torch.from_numpy(ds.labels.toarray().astype(np.float32))
    assert x.shape == (n, 2560) and y.shape == (n, len(terms))
    log('inputs loaded (dense, as the official loader holds them)', round(time.time() - t, 1), 's')

    forms = {}
    for name, path in axioms.items():
        nf1, nf2, nf3, nf4, rels, zero = load_normal_forms(path, terms)
        forms[name] = (tuple(torch.LongTensor(v).to(dev) for v in (nf1, nf2, nf3, nf4)), len(rels), len(zero))
        ops = 2 * len(nf1) + 3 * len(nf2) + 2 * len(nf3) + 2 * len(nf4)
        log(name, 'NF1-4', [len(nf1), len(nf2), len(nf3), len(nf4)], 'relations', len(rels), 'zero classes', len(zero),
            'class operands per step', ops)

    torch.manual_seed(0)
    loader = FastTensorDataLoader(x, y, batch_size=37, shuffle=True)
    t = time.time()
    it = iter(loader)  # the per-epoch shuffle: gather of the whole feature and label tensors
    shuffle_s = time.time() - t
    log('SHUFFLE seconds per epoch (CPU gather of features and labels)', round(shuffle_s, 2))
    batches = len(loader)

    def run(name, mode, count, net, optimizer, normal):
        """count batches in chunks; returns ms/batch per chunk."""
        out = []
        done = 0
        while done < count:
            torch.cuda.synchronize()
            t0 = time.time()
            k = min(CHUNK, count - done)
            for _ in range(k):
                bf, bl = next(it)
                bf, bl = bf.to(dev), bl.to(dev)
                optimizer.zero_grad()
                logits = net(bf)
                loss = torch.nn.functional.binary_cross_entropy(logits, bl)
                el_loss = net.el_loss(normal)
                total = loss + el_loss
                if mode == 'released':
                    loss.detach().item()
                    el_loss.detach().item()
                total.backward()
                optimizer.step()
            torch.cuda.synchronize()
            done += k
            out.append(1000 * (time.time() - t0) / k)
            log(f'CHUNK {name} {mode} batches {done}/{count} ms/batch {out[-1]:.1f}')
        return out

    results = []
    net = None
    order = ['paper', 'march', 'march', 'paper'] if 'march' in axioms else ['paper', 'paper']
    for name in order:
        normal, nrel, nzero = forms[name]
        torch.manual_seed(0)
        del net
        net = DeepGOModel(2560, len(terms), nzero, nrel, dev).to(dev)
        optimizer = torch.optim.Adam(net.parameters(), lr=5e-4)
        net.train()
        run(name, 'warmup', WARM, net, optimizer, normal)
        off = run(name, 'off', OFF, net, optimizer, normal)
        rel = run(name, 'released', REL, net, optimizer, normal)
        results.append((name, np.mean(off), np.mean(rel)))
        log(f'RESULT {name} off {np.mean(off):.2f} ms/batch (chunk range {min(off):.1f}-{max(off):.1f}), released {np.mean(rel):.2f}')

    for name in axioms:
        r = [v for v in results if v[0] == name]
        off, rel = np.mean([v[1] for v in r]), np.mean([v[2] for v in r])
        log(f'SUMMARY {name}: off {off:.2f} ms/batch, released {rel:.2f} ms/batch; epoch ({batches} batches) '
            f'off {off * batches / 1000 + shuffle_s:.1f} s, released {rel * batches / 1000 + shuffle_s:.1f} s (incl. shuffle {shuffle_s:.1f} s)')
    if 'march' in axioms:
        log('SUMMARY paper/march off',
            np.mean([v[1] for v in results if v[0] == 'paper']) / np.mean([v[1] for v in results if v[0] == 'march']),
            'released', np.mean([v[2] for v in results if v[0] == 'paper']) / np.mean([v[2] for v in results if v[0] == 'march']))

    # One development pass and the official per-epoch ROC AUC (validation_labels flattened, sklearn roc_curve on the CPU).
    dsd = load(root, 'development', sequences=False)
    dx = torch.from_numpy(np.ascontiguousarray(features.read(root, dsd), dtype=np.float32))
    dy = dsd.labels.toarray().astype(np.float32)
    net.eval()
    t = time.time()
    preds = []
    with torch.no_grad():
        for i in range(0, len(dx), 37):
            preds.append(net(dx[i:i + 37].to(dev)).cpu().numpy())
    preds = np.concatenate(preds)
    log('DEV PASS seconds', round(time.time() - t, 2), 'proteins', len(dx))
    t = time.time()
    roc = compute_roc(dy, preds)
    log('ROC AUC seconds per epoch (official compute_roc on the dev set)', round(time.time() - t, 2), 'auc', round(roc, 4))


if __name__ == '__main__':
    main()
