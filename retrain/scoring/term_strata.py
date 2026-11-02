"""The main gains on common and rare terms (Supplementary Note 10): test micro-AUPRC differences on the reviewed labels over
all terms, over the common terms (held by at least 1,000 HQ proteins, the line of the rare-function filter) and over the rest.
A stratum is the pairs of its terms; every protein and cluster is kept.

    python -m retrain.scoring.term_strata [--data DIR] [--device cuda|cpu]

A score is computed with the masked metrics' code (retrain.scoring.masked): the scorer's 10,000 resamples of the test clusters
give the paired 95% intervals, and a resample that draws every cluster once gives the point estimate. A model's value is the
mean over its seeds.

Writes results/term_strata/gains.tsv (stratum, terms, positives, method, training_set, minus_training_set, value, minus_value,
difference, ci_low, ci_high). Needs results/main. One minute on one GPU.
"""
import argparse

import numpy as np
import torch

from archipelago.evaluate import cluster_draws
from archipelago.go import Propagation, read_parents
from retrain.common import add_argument, data_root, load, log
from retrain.scoring.masked import replicate_auprc
from retrain.scoring.score import PUBLISHED, bundle, fits
from retrain.scoring.test_labels import annotation_sets

COMMON = 1000
MODELS = [('finch', s) for s in ('hq', 'hq_random', 'archipelago')] + [(m, s) for m in PUBLISHED for s in ('original', 'hq', 'archipelago')]
DIFFERENCES = [(m, 'archipelago', 'original') for m in PUBLISHED] + [(m, 'archipelago', 'hq') for m in ('finch', *PUBLISHED)] + \
              [('finch', 'archipelago', 'hq_random'), ('finch', 'hq_random', 'hq')]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--device', choices=['cuda', 'cpu'], help='default: cuda if available')
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    out = root / 'results' / 'term_strata'
    out.mkdir(parents=True, exist_ok=False)
    device = torch.device(a.device or ('cuda' if torch.cuda.is_available() else 'cpu'))
    torch.backends.cuda.matmul.allow_tf32 = False
    tensor = lambda x: torch.as_tensor(x, device=device)  # noqa: E731

    test, vocabulary, pairs, labels_of = annotation_sets(root)
    official = tensor(labels_of(pairs['all']))
    common = tensor(np.asarray(load(root, 'hq').labels.sum(axis=0)).ravel() >= COMMON)
    strata = {'all': torch.ones_like(common), 'common': common, 'rare': ~common}

    propagate = Propagation(read_parents(bundle(root) / 'go-basic.obo', vocabulary))
    cluster, draws = cluster_draws(test.columns['cluster'])
    cluster = tensor(cluster)
    draws = torch.as_tensor(draws.T, dtype=torch.float64, device=device)
    once = torch.ones((draws.shape[0], 1), dtype=torch.float64, device=device)  # every cluster once: the test set itself

    point, replicates = {}, {}
    for method, training_set, seeds in fits(root, 'main'):
        if (method, training_set) not in MODELS:
            continue
        per_seed = {}
        for files in seeds.values():
            s = propagate(torch.as_tensor(np.load(files['test']), dtype=torch.float32, device=device))
            for name, terms in strata.items():
                k = terms[None, :].expand_as(official)
                per_seed.setdefault(name, []).append((float(replicate_auprc(s, official, k, cluster, once)[0]),
                                                      replicate_auprc(s, official, k, cluster, draws)))
        for name, values in per_seed.items():
            point[name, method, training_set] = float(np.mean([v for v, _ in values]))
            replicates[name, method, training_set] = np.mean([r for _, r in values], 0)
        log(method, training_set)

    with open(out / 'gains.tsv', 'w', newline='\n') as handle:
        handle.write('stratum\tterms\tpositives\tmethod\ttraining_set\tminus_training_set\tvalue\tminus_value\tdifference\tci_low\t'
                     'ci_high\n')
        for name, terms in strata.items():
            for method, x, y in DIFFERENCES:
                low, high = np.percentile(replicates[name, method, x] - replicates[name, method, y], [2.5, 97.5])
                vx, vy = point[name, method, x], point[name, method, y]
                handle.write(f'{name}\t{int(terms.sum())}\t{int(official[:, terms].sum())}\t{method}\t{x}\t{y}\t{vx}\t{vy}\t'
                             f'{vx - vy}\t{low}\t{high}\n')


if __name__ == '__main__':
    main()
