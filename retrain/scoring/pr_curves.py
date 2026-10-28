"""Micro-averaged precision-recall curves on the test set (Fig. 1b and the Supplementary Information), from the scored
files of the main run.

    python -m retrain.scoring.pr_curves [--data DIR] [--device cuda|cpu]

The scorer's code path (archipelago.scorer): test labels and scores propagated with archipelago.go.Propagation, NaN
scores set to 0 and clipped to [0, 1]. All protein-term pairs are pooled and sorted by score; each distinct score is one
operating point (precision = TP / predicted, recall = TP / positives), so the step area is the scorer's micro-AUPRC. Each
curve is read at recall k/500, k = 1..500 (the precision of the first operating point reaching that recall; no
interpolation). Per model: one curve per seed, `mean` (the per-recall average of the seeds' curves; its exact area is the
seed mean of AUPRC) and, with several seeds, `ensemble` (the curve of the mean of the seeds' propagated scores).

Writes results/pr_curves/curves.tsv (method, training_set, seed, recall, precision) and check.tsv; stops unless every
exact step area equals the micro-AUPRC of results/main/metrics.tsv to three decimals.
"""
import argparse
import csv
import json

import numpy as np

from archipelago.go import Propagation, read_parents
from archipelago.io import read_labels, read_table, read_vocabulary
from retrain.common import add_argument, data_root, log

MODELS = {('finch', 'archipelago')} | {(m, t) for m in ('dpfunc', 'transfew', 'deepgose')
                                         for t in ('original', 'hq', 'archipelago')}
GRID = np.arange(1, 501) / 500


def main():
    import torch
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--device', choices=['cuda', 'cpu'])
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    out = root / 'results' / 'pr_curves'
    out.mkdir(parents=True, exist_ok=False)
    device = torch.device(a.device or ('cuda' if torch.cuda.is_available() else 'cpu'))
    config_path = root / 'scoring' / 'main.json'
    config = json.loads(config_path.read_text())
    data = config_path.parent / config['data']
    vocabulary = read_vocabulary(data / 'vocabulary.txt')
    ids = [r['protein_id'] for r in read_table(data / 'benchmark' / 'test.tsv')]
    propagate = Propagation(read_parents(data / 'go-basic.obo', vocabulary))
    y = propagate(read_labels(data / 'benchmark' / 'test_labels.tsv', ids, vocabulary).astype(np.float32)) > 0
    y = torch.as_tensor(y, device=device).reshape(-1)
    positives = int(y.sum())
    with open(root / 'results' / 'main' / 'metrics.tsv', newline='') as handle:
        scored = {(r['method'], r['training_set'], r['seed']): float(r['value'])
                  for r in csv.DictReader(handle, delimiter='\t') if r['split'] == 'test' and r['metric'] == 'auprc'}

    def load(path):  # as the scorer: rows in benchmark order, or reordered by the file's .ids.txt
        values = np.load(path)
        ids_file = path.with_suffix('.ids.txt')
        if ids_file.exists():
            position = {p: i for i, p in enumerate(ids_file.read_text().split())}
            values = values[[position[p] for p in ids]]
        if values.shape != (len(ids), len(vocabulary)):
            raise ValueError(f'{path}: shape {values.shape}')
        return torch.as_tensor(np.ascontiguousarray(values, dtype=np.float32), device=device)

    def curve(s):  # propagated scores -> (exact step area, precision at GRID)
        s = s.nan_to_num(0.0).clamp(0, 1).add(0.0).reshape(-1)
        s, order = torch.sort(s, descending=True, stable=True)
        tp = torch.cumsum(y[order].double(), 0)
        end = torch.nonzero(torch.cat([s[1:] != s[:-1], torch.tensor([True], device=device)])).squeeze(1)
        tp = tp[end]
        precision, recall = tp / (end + 1).double(), tp / positives
        area = float((torch.diff(tp, prepend=tp.new_zeros(1)) / positives * precision).sum())
        at = torch.searchsorted(recall, torch.as_tensor(GRID, device=device) - 1e-12)
        return area, precision[at].cpu().numpy()

    rows, checks = [], []
    for model in config['models']:
        key = (model['method'], model['training_set'])
        if key not in MODELS:
            continue
        curves, areas, total = [], [], None
        for seed, files in model['seeds'].items():
            s = propagate(load(config_path.parent / files['test']), inplace=True)
            area, precision = curve(s)
            curves.append((seed, area, precision))
            areas.append(area)
            total = s.double() if total is None else total.add_(s)
        curves.append(('mean', float(np.mean(areas)), np.mean([c[2] for c in curves], 0)))
        if len(areas) > 1:
            curves.append(('ensemble',) + curve(total.div_(len(areas)).float()))
        for seed, area, precision in curves:
            rows += [(*key, seed, f'{r:.4f}', f'{p:.6f}') for r, p in zip(GRID, precision)]
            expected = scored[(*key, seed)]
            checks.append((*key, seed, f'{area:.6f}', f'{expected:.6f}', round(area, 3) == round(expected, 3)))
            log(*checks[-1])
    with open(out / 'curves.tsv', 'w', newline='\n') as handle:
        handle.write('method\ttraining_set\tseed\trecall\tprecision\n')
        handle.writelines('\t'.join(r) + '\n' for r in rows)
    with open(out / 'check.tsv', 'w', newline='\n') as handle:
        handle.write('method\ttraining_set\tseed\tarea\tscored_auprc\tsame_3_decimals\n')
        handle.writelines('\t'.join(map(str, c)) + '\n' for c in checks)
    if not all(c[-1] for c in checks):
        raise SystemExit('a curve area differs from the scored micro-AUPRC (see check.tsv)')


if __name__ == '__main__':
    main()
