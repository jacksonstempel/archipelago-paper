"""The main gains on common and rare terms (Supplementary Note 10): test micro-AUPRC differences on the reviewed labels over
all terms, over the common terms (held by at least 1,000 HQ proteins, the line of the rare-function filter) and over the rest.
A stratum is the pairs of its terms; every protein and cluster is kept.

    python -m retrain.scoring.term_strata [--data DIR] [--labels reviewed|masked] [--device cuda|cpu]

A score is an exact micro-AUPRC over the pairs of a stratum (replicate_auprc): the scorer's 10,000 resamples of the test
clusters give the paired 95% intervals, and a resample that draws every cluster once gives the point estimate. A model's value
is the mean over its seeds.

--labels masked: the masked test labels instead (the reviewed GOA assertions are neither positive nor negative, in no stratum;
the strata are the same terms), written to results_masked/term_strata/gains.tsv.

Writes results/term_strata/gains.tsv (stratum, terms, positives, method, training_set, minus_training_set, value, minus_value,
difference, ci_low, ci_high). Needs results/main. One minute on one GPU.
"""
import argparse

import numpy as np
import torch

from archipelago.evaluate import cluster_draws
from archipelago.go import Propagation, read_parents
from retrain.common import add_argument, data_root, load, log
from retrain.scoring.score import PUBLISHED, RESULTS, fits, results_dir
from retrain.scoring.test_labels import annotation_sets, version_labels

COMMON = 1000
MODELS = [('finch', s) for s in ('hq', 'hq_random', 'archipelago')] + [(m, s) for m in PUBLISHED for s in ('original', 'hq', 'archipelago')]
DIFFERENCES = [(m, 'archipelago', 'original') for m in PUBLISHED] + [(m, 'archipelago', 'hq') for m in ('finch', *PUBLISHED)] + \
              [('finch', 'archipelago', 'hq_random'), ('finch', 'hq_random', 'hq')]


def replicate_auprc(s, y, keep, cluster, draws, chunk=1000):
    """Micro-AUPRC of propagated scores s against labels y over the pairs in keep in every bootstrap replicate. cluster:
    the cluster index of each protein (tensor); draws: (clusters x replicates, float64) the times each cluster was drawn.
    Each positive takes the weighted precision at the end of its tied-score group: the pairs are counted per cluster in
    segments that end at the positives' group ends, so the weighted sums are (segments x clusters) @ (clusters x
    replicates)."""
    device, terms = s.device, s.shape[1]
    flat = torch.nonzero(keep.reshape(-1)).squeeze(1)
    order = torch.argsort(s.reshape(-1)[flat], descending=True, stable=True)
    flat = flat[order]
    sv, yv = s.reshape(-1)[flat], y.reshape(-1)[flat]
    ends = torch.nonzero(torch.cat([sv[1:] != sv[:-1], torch.ones(1, dtype=torch.bool, device=device)])).squeeze(1)
    group_end = ends[torch.searchsorted(ends, torch.arange(len(sv), device=device))]
    segment_ends = torch.unique(group_end[yv])  # sorted
    last = int(segment_ends[-1]) + 1  # pairs after the last positive's group do not enter the AUPRC
    segment = torch.searchsorted(segment_ends, torch.arange(last, device=device))
    of = cluster[flat[:last] // terms]
    k, segments = draws.shape[0], len(segment_ends)
    every = torch.bincount(segment * k + of, minlength=segments * k).reshape(segments, k).double()
    at = yv[:last]
    positive = torch.sparse_coo_tensor(torch.stack([segment[at], of[at]]),
                                       torch.ones(int(at.sum()), dtype=torch.float64, device=device),
                                       (segments, k), check_invariants=True).coalesce()  # repeated entries summed
    out = []
    for a in range(0, draws.shape[1], chunk):
        w = draws[:, a:a + chunk]
        p = torch.sparse.mm(positive, w)
        tp, n = torch.cumsum(p, 0), torch.cumsum(every @ w, 0)
        out.append(torch.where(n > 0, p * tp / n.clamp_min(1e-12), torch.zeros_like(p)).sum(0) / p.sum(0))
    return torch.cat(out).cpu().numpy()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--device', choices=['cuda', 'cpu'], help='default: cuda if available')
    parser.add_argument('--labels', choices=list(RESULTS), default='reviewed', help='the version of the main test labels')
    add_argument(parser)
    a = parser.parse_args(argv)
    root = data_root(a.data)
    out = results_dir(root, a.labels) / 'term_strata'
    out.mkdir(parents=True, exist_ok=False)
    device = torch.device(a.device or ('cuda' if torch.cuda.is_available() else 'cpu'))
    torch.backends.cuda.matmul.allow_tf32 = False
    tensor = lambda x: torch.as_tensor(x, device=device)  # noqa: E731

    test, vocabulary, pairs, labels_of = annotation_sets(root)
    if a.labels == 'masked':
        official, left_out = (tensor(m) for m in version_labels(pairs, labels_of, 'masked'))
    else:
        official, left_out = tensor(labels_of(pairs['all'])), tensor(np.zeros(test.labels.shape, dtype=bool))
    counted = ~left_out  # the pairs that count
    common = tensor(np.asarray(load(root, 'hq').labels.sum(axis=0)).ravel() >= COMMON)
    strata = {'all': torch.ones_like(common), 'common': common, 'rare': ~common}

    propagate = Propagation(read_parents(root / 'datasets' / 'ontology' / 'go-basic.obo', vocabulary))
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
                k = counted & terms[None, :]
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
