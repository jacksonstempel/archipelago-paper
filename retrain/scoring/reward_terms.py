"""The main gains on the terms that the split search rewarded for computational support, and on all other terms
(Supplementary Note 10 and its table), on the main test labels.

    python -m retrain.scoring.reward_terms [--data DIR] [--device cuda|cpu]

The benchmark split was searched with a reward of which a fifth went to GO terms that a computational training candidate
(a Swiss-Prot protein with an IBA annotation or a TrEMBL protein with an IEA annotation) carries beyond its high-quality
labels and that keep such a protein in training (Supplementary Note 10). rewarded_terms.tsv, beside this file, lists every
term with an opportunity weight (column `rewarded` = 1: the term made a nonzero contribution to that part of the reward in
the chosen split; `retained_extra_proteins` > 0: it has the weight and keeps computational support in training). It was
recomputed from the split search's saved state, which reproduced the search's own record of that reward exactly; the
search and its state are not part of the release.

Strata (scoring vocabulary, 8,739 terms; a stratum is the pairs of its terms, as in retrain.scoring.term_strata; all 5,000
test proteins and their clusters are kept):
  rewarded / other           the 2,327 rewarded terms in the vocabulary / the other 6,412
  eligible / other_eligible  sensitivity: every term with the weight and retained computational support / the rest
  all                        check: reproduces results_masked/term_strata/gains.tsv to 1e-9

A score is the exact micro-AUPRC over a stratum's pairs on the main test labels (the pairs that the reviewed GOA assertions
add are left out of every stratum), with the scorer's 10,000 cluster-bootstrap resamples for the paired intervals
(retrain.scoring.term_strata.replicate_auprc); a model's value and replicates are means over its runs. Models: the published
methods on original, hq, hq_explabels and archipelago and Finch on hq, hq_random and archipelago, from scoring/main.json and
scoring/explabels_training.json (retrain.scoring.score; a model in both must have the same test files).

Writes results_masked/reward_terms/gains.tsv (stratum, terms, positives, method, training_set, minus_method,
minus_training_set, value, minus_value, difference, ci_low, ci_high). Needs results_masked/term_strata (retrain.scoring.
term_strata --labels masked). About two minutes on one GPU.
"""
import argparse
import csv
import json
import os
import time
from pathlib import Path

import numpy as np
import torch

from archipelago.evaluate import cluster_draws
from archipelago.go import Propagation, read_parents
from retrain.common import add_argument, data_root, load, log
from retrain.scoring.score import bundle, results_dir
from retrain.scoring.term_strata import replicate_auprc
from retrain.scoring.test_labels import annotation_sets, version_labels

TERMS = Path(__file__).with_name('rewarded_terms.tsv')
PUBLISHED = ('transfew', 'deepgose', 'dpfunc')
MODELS = [(m, s) for m in PUBLISHED for s in ('original', 'hq', 'hq_explabels', 'archipelago')] + \
         [('finch', s) for s in ('hq', 'hq_random', 'archipelago')]
PAIRS = [((m, a), (m, b)) for m in PUBLISHED for a, b in (('archipelago', 'original'), ('hq', 'original'),
                                                          ('hq', 'hq_explabels'), ('archipelago', 'hq'))] + \
        [(('finch', 'archipelago'), ('finch', 'hq_random')), (('finch', 'archipelago'), ('finch', 'hq'))] + \
        [(('finch', 'archipelago'), (m, 'archipelago')) for m in PUBLISHED]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--device', choices=['cuda', 'cpu'], help='default: cuda if available')
    add_argument(parser)
    a = parser.parse_args(argv)
    root = data_root(a.data)
    out = results_dir(root, 'masked') / 'reward_terms'
    reference_path = results_dir(root, 'masked') / 'term_strata' / 'gains.tsv'
    if not reference_path.exists():
        raise FileNotFoundError(f'{reference_path}: run retrain.scoring.term_strata --labels masked first')
    out.mkdir(parents=True, exist_ok=False)
    device = torch.device(a.device or ('cuda' if torch.cuda.is_available() else 'cpu'))
    torch.backends.cuda.matmul.allow_tf32 = False
    tensor = lambda x: torch.as_tensor(x, device=device)  # noqa: E731

    # ---- the main test labels and the models' files
    test, vocabulary, pairs, labels_of = annotation_sets(root)
    y, left_out = (tensor(m) for m in version_labels(pairs, labels_of, 'masked'))
    keep = ~left_out
    files = {}
    for run in ('main', 'explabels_training'):
        config = root / 'scoring' / f'{run}.json'
        for m in json.loads(config.read_text())['models']:
            key = (m['method'], m['training_set'])
            tests = {seed: os.path.realpath(config.parent / f['test']) for seed, f in m['seeds'].items()}
            if files.setdefault(key, tests) != tests:
                raise ValueError(f'{key} has different test files in {config}')
    missing = [k for k in MODELS if k not in files]
    if missing:
        raise ValueError(f'models missing from the configs: {missing}')

    # ---- the strata
    rows = list(csv.DictReader(open(TERMS, newline=''), delimiter='\t'))
    index = {t: j for j, t in enumerate(vocabulary)}

    def term_mask(names):
        v = torch.zeros(len(vocabulary), dtype=torch.bool, device=device)
        v[[index[t] for t in names if t in index]] = True
        return v

    rewarded = term_mask([r['go_id'] for r in rows if r['rewarded'] == '1'])
    eligible = term_mask([r['go_id'] for r in rows if int(r['retained_extra_proteins']) > 0])
    common = tensor(np.asarray(load(root, 'hq').labels.sum(axis=0)).ravel() >= 1000)
    strata = {'all': torch.ones_like(rewarded), 'rewarded': rewarded, 'other': ~rewarded,
              'eligible': eligible, 'other_eligible': ~eligible}
    pos = lambda t: int((y & keep)[:, t].sum())  # noqa: E731
    log({k: (int(v.sum()), pos(v)) for k, v in strata.items()}, '(terms, positive pairs, main labels)')
    log(f'common terms {int(common.sum())}: rewarded and common {int((rewarded & common).sum())} '
        f'(positives {pos(rewarded & common)}); eligible and common {int((eligible & common).sum())}; '
        f'rewarded terms with a test positive {int(((y & keep)[:, rewarded]).any(0).sum())}')

    # ---- every model on every stratum
    propagate = Propagation(read_parents(bundle(root) / 'go-basic.obo', vocabulary))
    cluster, draws = cluster_draws(test.columns['cluster'])
    cluster = tensor(cluster)
    draws = torch.as_tensor(draws.T, dtype=torch.float64, device=device)
    once = torch.ones((draws.shape[0], 1), dtype=torch.float64, device=device)  # every cluster once: the test set itself
    point, replicates, start = {}, {}, time.time()
    for key in MODELS:
        per_seed = {}
        for path in files[key].values():
            s = propagate(torch.as_tensor(np.load(path), dtype=torch.float32, device=device))
            for name, terms in strata.items():
                k = keep & terms[None, :]
                per_seed.setdefault(name, []).append((float(replicate_auprc(s, y, k, cluster, once)[0]),
                                                      replicate_auprc(s, y, k, cluster, draws)))
        for name, values in per_seed.items():
            point[(name, *key)] = float(np.mean([v for v, _ in values]))
            replicates[(name, *key)] = np.mean([r for _, r in values], 0)
        log(f'{time.time() - start:6.0f} s', *key, len(per_seed['all']), 'runs',
            {n: round(point[(n, *key)], 4) for n in strata})

    path = out / 'gains.tsv'
    with open(path, 'w', newline='\n') as handle:
        handle.write('stratum\tterms\tpositives\tmethod\ttraining_set\tminus_method\tminus_training_set\tvalue\tminus_value\t'
                     'difference\tci_low\tci_high\n')
        for name, terms in strata.items():
            for x, z in PAIRS:
                vx, vz = point[(name, *x)], point[(name, *z)]
                low, high = np.percentile(replicates[(name, *x)] - replicates[(name, *z)], [2.5, 97.5])
                handle.write(f'{name}\t{int(terms.sum())}\t{pos(terms)}\t{x[0]}\t{x[1]}\t{z[0]}\t{z[1]}\t{vx!r}\t{vz!r}\t'
                             f'{vx - vz!r}\t{low!r}\t{high!r}\n')

    # ---- check: the all-terms rows reproduce the term strata of the main labels
    mine = {(r['method'], r['training_set'], r['minus_training_set']): r for r in csv.DictReader(open(path), delimiter='\t')
            if r['stratum'] == 'all' and r['method'] == r['minus_method']}
    checked = 0
    for r in csv.DictReader(open(reference_path), delimiter='\t'):
        key = (r['method'], r['training_set'], r['minus_training_set'])
        if r['stratum'] == 'all' and key in mine:
            for c in ('value', 'minus_value', 'difference', 'ci_low', 'ci_high'):
                if abs(float(mine[key][c]) - float(r[c])) > 1e-9:
                    raise ValueError(f'{key} {c}: {mine[key][c]} vs {reference_path} {r[c]}')
            checked += 1
    log(f'all-terms rows: {checked} pairs reproduce {reference_path}')
    if checked < 5:
        raise ValueError('too few reference pairs checked')
    log('wrote', path)


if __name__ == '__main__':
    main()
