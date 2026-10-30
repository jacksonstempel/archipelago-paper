"""Scoring with the reviewed GOA assertions masked (main-text Table 2, Supplementary Information). Positives are the test
labels without any reviewed assertion (the experimental and curated labels the benchmark was made with); the pairs that
the reviewed assertions add, accepted or not, with their parent terms, count as neither positive nor negative. Every test
protein is kept.

    python -m retrain.scoring.masked [--data DIR] [--device cuda|cpu]

The scorer cannot leave pairs out (and its propagation would undo a zeroed score), so the metrics are computed here with
the scorer's definitions (archipelago.evaluate) over the remaining pairs. Every fit of the main run: raw test scores
propagated as the scorer propagates them; micro-AUPRC (exact average precision, tied scores grouped), micro-Fmax
(thresholds 0.01 to 1.00) and Smin (0.00 to 1.00; information accretion from the HQ labels of the scoring bundle, as the
scorer estimates it); seed means. Paired 95% intervals of micro-AUPRC for the models of the paper's label-set tables
(BOOTSTRAP), with the scorer's 10,000 resamples of the test clusters (archipelago.evaluate.cluster_draws, its fixed seed):
a protein weighs as often as its cluster was drawn, and a model's replicate is the mean over its seeds.

Checks, with the official labels and no mask (against results/main): the first Finch fit reproduces the scorer's test
AUPRC, and Finch on ArchipelaGO minus Finch on HQ + random reproduces the scorer's interval.

Writes results/masked_test_labels/: metrics.tsv (method, training_set, seed, split, metric, value; split `test`, every seed
and their `mean`), contrasts.tsv (micro-AUPRC, every pair of BOOTSTRAP models with the same method or training set;
columns as retrain.scoring.score) and counts.tsv (proteins scored, proteins with a positive, clusters, positive and masked
pairs, terms with a positive, for the official labels, `all`, and the masked ones, `masked`). Needs results/main. A few
minutes on one GPU.
"""
import argparse
import itertools

import numpy as np
import torch

from archipelago.evaluate import FMAX_THRESHOLDS, SMIN_THRESHOLDS, cluster_draws, paired_difference
from archipelago.go import Propagation, information_accretion, read_parents
from archipelago.io import read_labels, read_table
from retrain.common import add_argument, data_root, log
from retrain.scoring.score import PUBLISHED, bundle, fits
from retrain.scoring.test_labels import annotation_sets, write_counts

BOOTSTRAP = [(m, s) for m in PUBLISHED for s in ('original', 'goa', 'hq', 'archipelago')] + \
            [('finch', s) for s in ('goa', 'hq', 'hq_random', 'archipelago')]
METRICS = ('auprc', 'fmax', 'smin')


def metrics(s, y, keep, ia):
    """micro-AUPRC, micro-Fmax and Smin of propagated scores s (proteins x terms, float32 tensor) against labels y over the
    pairs in keep (boolean tensors of the same shape; no positive may be left out), as archipelago.evaluate defines them;
    ia: information accretion of each term (float64 tensor). At threshold t a pair is predicted when its score is >= t (t
    rounded to float32); Smin's misinformation counts only the kept pairs."""
    if (y & ~keep).any():
        raise ValueError('a positive pair is masked')
    flat = torch.nonzero(keep.reshape(-1)).squeeze(1)
    order = torch.argsort(s.reshape(-1)[flat], descending=True, stable=True)
    flat = flat[order]
    sv, yv = s.reshape(-1)[flat], y.reshape(-1)[flat]
    tp = torch.cumsum(yv, 0, dtype=torch.float64)
    ends = torch.nonzero(torch.cat([sv[1:] != sv[:-1], torch.ones(1, dtype=torch.bool, device=sv.device)])).squeeze(1)
    group_end = ends[torch.searchsorted(ends, torch.arange(len(sv), device=sv.device))]
    positives = tp[-1]
    auprc = float((tp[group_end] / (group_end + 1.0))[yv].sum() / positives)

    ascending = sv.flip(0)

    def first(cumulative, thresholds):
        """The cumulative sum over the pairs predicted at each threshold (the first ones in score order)."""
        t = torch.as_tensor(np.asarray(thresholds, dtype=np.float32), device=sv.device)
        n = len(sv) - torch.searchsorted(ascending, t, side='left')
        return torch.where(n > 0, cumulative[(n - 1).clamp_min(0)], torch.zeros_like(cumulative[:1])), n

    hits, predicted = first(tp, FMAX_THRESHOLDS)
    fmax = float((2 * hits / (predicted + positives)).max())
    weight = ia[flat % s.shape[1]]
    found, _ = first(torch.cumsum(weight * yv, 0), SMIN_THRESHOLDS)
    wrong, _ = first(torch.cumsum(weight * ~yv, 0), SMIN_THRESHOLDS)
    remaining = (weight * yv).sum() - found
    smin = float(torch.hypot(remaining / s.shape[0], wrong / s.shape[0]).min())
    return dict(auprc=auprc, fmax=fmax, smin=smin)


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


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--device', choices=['cuda', 'cpu'], help='default: cuda if available')
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    out = root / 'results' / 'masked_test_labels'
    if out.exists():
        raise FileExistsError(f'{out} exists')
    device = torch.device(a.device or ('cuda' if torch.cuda.is_available() else 'cpu'))
    torch.backends.cuda.matmul.allow_tf32 = False

    test, vocabulary, pairs, labels_of = annotation_sets(root)
    official, positive = labels_of(pairs['all']), labels_of(pairs['other'])
    masked = labels_of(pairs['all_goa']) & ~positive
    counts = []
    for name, y, mask in (('all', official, np.zeros_like(official)), ('masked', positive, masked)):
        counts.append(dict(variant=name, proteins=len(y), proteins_with_positive=int(y.any(1).sum()),
                           clusters=len(set(test.columns['cluster'])), propagated_pairs=int(y.sum()),
                           masked_pairs=int(mask.sum()), terms_with_positive=int(y.any(0).sum())))
        log(counts[-1])

    data = bundle(root)
    parents = read_parents(data / 'go-basic.obo', vocabulary)
    propagate = Propagation(parents)
    hq_ids = [r['protein_id'] for r in read_table(data / 'archipelago' / 'proteins.tsv') if r['part'] == 'hq']
    hq = propagate(read_labels(data / 'archipelago' / 'labels.tsv', hq_ids, vocabulary, other_proteins=True)
                   .astype(np.float32)) > 0
    ia = torch.as_tensor(information_accretion(hq, parents), dtype=torch.float64, device=device)
    cluster, draws = cluster_draws(test.columns['cluster'])
    cluster = torch.as_tensor(cluster, device=device)
    draws = torch.as_tensor(draws.T, dtype=torch.float64, device=device)

    def scores(path):
        s = np.load(path)
        if s.shape != official.shape:
            raise ValueError(f'{path}: scores of shape {s.shape}, expected {official.shape}')
        return propagate(torch.as_tensor(s, dtype=torch.float32, device=device))

    models = fits(root, 'main')
    seeds_of = {(m, t): seeds for m, t, seeds in models}

    # checks: the official labels without a mask reproduce the scorer (results/main)
    y, keep = torch.as_tensor(official, device=device), torch.ones(official.shape, dtype=torch.bool, device=device)
    main_results = root / 'results' / 'main'
    scorer = {(r['method'], r['training_set'], r['seed']): float(r['value']) for r in read_table(main_results / 'metrics.tsv')
              if r['split'] == 'test' and r['metric'] == 'auprc'}
    method, training_set, seeds = next(m for m in models if m[0] == 'finch')
    seed, files = next(iter(seeds.items()))
    found, wanted = metrics(scores(files['test']), y, keep, ia)['auprc'], scorer[method, training_set, seed]
    log(f'check {method}/{training_set}/seed{seed}: test AUPRC {found:.10f}, scorer {wanted:.10f}')
    if abs(found - wanted) > 1e-6:
        raise ValueError('the official labels without a mask do not reproduce the scorer')
    finch = {t: np.mean([replicate_auprc(scores(f['test']), y, keep, cluster, draws)
                         for f in seeds_of['finch', t].values()], 0) for t in ('archipelago', 'hq_random')}
    low, high = np.percentile(finch['archipelago'] - finch['hq_random'], [2.5, 97.5])
    row = next(r for r in read_table(main_results / 'contrasts.tsv') if r['metric'] == 'auprc' and (
        r['method'], r['training_set'], r['minus_method'], r['minus_training_set']) == ('finch', 'archipelago', 'finch',
                                                                                      'hq_random'))
    log(f"check interval of finch archipelago - hq_random: [{low:.6f}, {high:.6f}], scorer [{float(row['ci_low']):.6f}, "
        f"{float(row['ci_high']):.6f}]")
    if max(abs(low - float(row['ci_low'])), abs(high - float(row['ci_high']))) > 1e-5:
        raise ValueError('the bootstrap does not reproduce the scorer interval')

    # every fit with the mask
    y, keep = torch.as_tensor(positive, device=device), torch.as_tensor(~masked, device=device)
    point, replicates, rows = {}, {}, []
    for method, training_set, seeds in models:
        per_seed, per_seed_replicates = [], []
        for seed, files in seeds.items():
            s = scores(files['test'])
            per_seed.append(metrics(s, y, keep, ia))
            if (method, training_set) in BOOTSTRAP:
                per_seed_replicates.append(replicate_auprc(s, y, keep, cluster, draws))
            rows += [(method, training_set, seed, k, v) for k, v in per_seed[-1].items()]
        mean = {k: float(np.mean([r[k] for r in per_seed])) for k in METRICS}
        rows += [(method, training_set, 'mean', k, v) for k, v in mean.items()]
        point[method, training_set] = mean['auprc']
        if per_seed_replicates:
            replicates[method, training_set] = np.mean(per_seed_replicates, 0)
        log(method, training_set, {k: round(v, 4) for k, v in mean.items()})
    out.mkdir(parents=True)
    with open(out / 'metrics.tsv', 'w', newline='\n') as handle:
        handle.write('method\ttraining_set\tseed\tsplit\tmetric\tvalue\n')
        handle.writelines(f'{m}\t{t}\t{seed}\ttest\t{k}\t{v!r}\n' for m, t, seed, k, v in rows)
    with open(out / 'contrasts.tsv', 'w', newline='\n') as handle:
        handle.write('metric\tmethod\ttraining_set\tminus_method\tminus_training_set\tdifference\tci_low\tci_high\n')
        for a, b in itertools.permutations(BOOTSTRAP, 2):
            if a[0] == b[0] or a[1] == b[1]:
                difference, (low, high) = paired_difference((point[a], replicates[a]), (point[b], replicates[b]))
                handle.write(f'auprc\t{a[0]}\t{a[1]}\t{b[0]}\t{b[1]}\t{difference!r}\t{float(low)!r}\t{float(high)!r}\n')
    write_counts(out / 'counts.tsv', counts)
    log(f'wrote {out}')


if __name__ == '__main__':
    main()
