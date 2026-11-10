"""Direct-evidence test labels and every model on every version of the test labels (Supplementary Note 2, main-text Table 2).

    python -m retrain.scoring.label_grid [--data DIR] [--device cuda|cpu]

The versions, built from the test annotations and the label review (retrain.scoring.test_labels.annotation_sets; the
review's release/benchmark_reviews/test.parquet says for each accepted assertion whether the review found experimental
support for it, support_type `experimental`, or accepted it as a supported inference, `inference`), all on the 5,000 test
proteins unless stated:

  experimental      the test annotations with an experimental evidence code; proteins without one are left out
  direct_negatives  positives: the test annotations other than the review's plus the assertions accepted with experimental
                    support; the official positives that only the assertions accepted as a supported inference give (with
                    their parent terms) count as negatives
  direct_masked     direct evidence only: the same positives; those pairs count as neither positive nor negative. Every
                    other pair is labelled as in the official (reviewed) labels
  masked            positives: the test annotations other than the review's; the pairs that every assertion the review
                    considered adds, accepted or not, count as neither positive nor negative (the main test labels)
  masked_negatives  computational annotations negative: the positives of masked; those pairs count as negatives (all 5,000
                    proteins; main-text Table 2, not a column of the grid)
  reviewed          the official labels
  all_goa           the test annotations other than the review's plus every assertion the review considered
  no_hit_0.001, no_hit_10   the official labels (in results_masked/label_grid: the masked ones) on the test proteins
                    without a training hit at E <= 1e-3 and E <= 10 (the paper's subsets, release/benchmark_no_hit)

A score is the exact micro-AUPRC over the pairs a version keeps (retrain.scoring.term_strata.replicate_auprc: a resample
that draws every cluster once gives the point estimate; the scorer's 10,000 resamples of the test clusters give paired
intervals), of the test scores propagated over the ontology; a model's value is the mean over its runs (DeepGO-SE: its
ensemble). Models: scoring/main.json and scoring/explabels_training.json (retrain.scoring.score; a model in both must have
the same test files). Writes:

1. results/direct_masked_test_labels/, results/direct_negatives_test_labels/ and results/masked_negatives_test_labels/:
   models.tsv (the published methods on original, goa, hq and archipelago and Finch on goa, hq, hq_random and archipelago;
   micro-AUPRC per run and mean) and contrasts.tsv (every pair of them within a method or a training set: the difference of
   the means and its paired 95% interval; main-text Table 2 reads the direct_masked and masked_negatives ones).
2. results/label_grid/grid.tsv (the published methods on original, hq_explabels, goa, hq and archipelago and Finch on the
   three original sets, hq_explabels, goa, hq, hq_random and archipelago; one column per version), versions.tsv (proteins,
   clusters, positive and left-out pairs, terms with a positive) and check.tsv: every grid value that a scored results
   folder has for the same version and model (main, explabels_training, experimental_test_labels, masked_test_labels,
   all_goa_test_labels, no_hit/E0.001, no_hit/E10, and the direct-evidence folders of 1.) must agree to 1e-6, or the run
   stops.
3. results_masked/label_grid/: the same with the two no-hit columns (and their rows of versions.tsv) on the main test
   labels instead of the reviewed ones, as in Supplementary Table 6 and Supplementary Figure 1; those columns are checked
   against results_masked/no_hit/E0.001 and E10.

Checks first: with the official labels, the first Finch run reproduces its micro-AUPRC in results/main, and Finch
ArchipelaGO - HQ + random its paired interval. Needs the results of retrain.scoring.score main and explabels_training,
retrain.scoring.test_labels experimental, masked and all_goa, retrain.scoring.no_hit, and retrain.scoring.headline
(results_masked/no_hit). A few minutes on one GPU.
"""
import argparse
import collections
import csv
import json
import os
import time
from pathlib import Path

import numpy as np
import torch

from archipelago.evaluate import cluster_draws
from archipelago.go import Propagation, read_parents
from retrain.common import add_argument, data_root, log
from retrain.scoring.no_hit import released_subsets
from retrain.scoring.score import ORIGINAL, results_dir
from retrain.scoring.term_strata import replicate_auprc
from retrain.scoring.test_labels import annotation_sets

PUBLISHED = ('deepgose', 'transfew', 'dpfunc')
TABLES = ('direct_masked', 'direct_negatives', 'masked_negatives')
BOOTSTRAP = [(m, s) for m in ('transfew', 'deepgose', 'dpfunc') for s in ('original', 'goa', 'hq', 'archipelago')] + \
            [('finch', s) for s in ('goa', 'hq', 'hq_random', 'archipelago')]
GRID_ROWS = [(m, s) for m in PUBLISHED for s in ('original', 'hq_explabels', 'goa', 'hq', 'archipelago')] + \
            [('finch', s) for s in (ORIGINAL['deepgose'], ORIGINAL['transfew'], ORIGINAL['dpfunc'], 'hq_explabels', 'goa', 'hq',
                                    'hq_random', 'archipelago')]
GRID_VERSIONS = ('experimental', 'direct_negatives', 'direct_masked', 'masked', 'reviewed', 'all_goa', 'no_hit_0.001',
                 'no_hit_10')
REFERENCES = (('reviewed', 'main'), ('reviewed', 'explabels_training'), ('experimental', 'experimental_test_labels'),
              ('masked', 'masked_test_labels'), ('all_goa', 'all_goa_test_labels'), ('no_hit_0.001', 'no_hit/E0.001'),
              ('no_hit_10', 'no_hit/E10'), ('direct_masked', 'direct_masked_test_labels'),
              ('direct_negatives', 'direct_negatives_test_labels'))  # (version, results folder) of the grid's check


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--device', choices=['cuda', 'cpu'], help='default: cuda if available')
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    results, results_masked = results_dir(root, 'reviewed'), results_dir(root, 'masked')
    outs = {name: results / f'{name}_test_labels' for name in TABLES}
    outs['grid'], outs['main_grid'] = results / 'label_grid', results_masked / 'label_grid'
    for folder in outs.values():
        if folder.exists():
            raise FileExistsError(f'{folder} exists')
    dev = torch.device(a.device or ('cuda' if torch.cuda.is_available() else 'cpu'))
    torch.backends.cuda.matmul.allow_tf32 = False

    # ---- labels
    test, vocabulary, pairs, closure = annotation_sets(root)
    ids, row_of = test.seq_ids, {h: i for i, h in enumerate(test.seq_ids)}
    official, base = closure(pairs['all']), closure(pairs['other'])
    masked = closure(pairs['all_goa']) & ~base
    log(f'positives {int(base.sum())} on {int(base.any(1).sum())} proteins; masked pairs {int(masked.sum())}; '
        f'official positives {int(official.sum())}')
    direct = closure(pairs['direct'])
    inferred = official & ~direct  # the official positives that only the inference-based acceptances give
    if (direct & ~official).any():
        raise ValueError('the direct-evidence labels have a positive that the official labels lack')
    inference = sorted(set(pairs['all']) - set(pairs['direct']))  # the assertions accepted as a supported inference
    by_inference = closure(pairs['other'] + inference) & ~base
    log(f"direct evidence: positives {int(direct.sum())} on {int(direct.any(1).sum())} proteins; the inference-based "
        f"acceptances add {int(by_inference.sum())} propagated pairs to the other rows, {int((by_inference & direct).sum())} "
        f"of them also given by the other rows or the experimental acceptances; {int(inferred.sum())} official positives on "
        f"{int(inferred.any(1).sum())} proteins only theirs (left out / negative)")

    propagate = Propagation(read_parents(root / 'datasets' / 'ontology' / 'go-basic.obo', vocabulary))
    cluster, draws = cluster_draws(test.columns['cluster'])
    cluster_of = torch.as_tensor(cluster, device=dev)
    draws = torch.as_tensor(draws.T, dtype=torch.float64, device=dev)  # clusters x replicates
    once = torch.ones((draws.shape[0], 1), dtype=torch.float64, device=dev)  # every cluster once: the test set itself

    def tensor(x):
        return torch.as_tensor(x, device=dev)

    def scores_of(path):
        values = np.load(path).astype(np.float32)
        if values.shape != (len(ids), len(vocabulary)):
            raise ValueError(f'{path}: shape {values.shape}')
        return propagate(torch.as_tensor(values, device=dev))

    def auprc(s, y, keep):
        return float(replicate_auprc(s, y, keep, cluster_of, once)[0])

    # ---- models (a model in both configs once, if its test files are the same)
    files = {}
    for run in ('main', 'explabels_training'):
        config = root / 'scoring' / f'{run}.json'
        c = json.loads(config.read_text())
        for m in c['models']:
            key = (m['method'], m['training_set'])
            tests = {seed: os.path.realpath(config.parent / f['test']) for seed, f in m['seeds'].items()}
            if files.setdefault(key, tests) != tests:
                raise ValueError(f'{key} has different test files in {config}')
    missing = [k for k in GRID_ROWS + BOOTSTRAP if k not in files]
    if missing:
        raise ValueError(f'models missing from the configs: {missing}')

    # ---- checks: the official labels reproduce the scorer's results/main
    y_official, keep_all = tensor(official), torch.ones(official.shape, dtype=torch.bool, device=dev)
    reference = {(r['method'], r['training_set'], r['seed']): r for r in csv.DictReader(
        open(results / 'main' / 'metrics.tsv', newline=''), delimiter='\t') if r['split'] == 'test' and r['metric'] == 'auprc'}
    seed, path = next(iter(files['finch', 'goa'].items()))
    value, want = auprc(scores_of(path), y_official, keep_all), float(reference['finch', 'goa', seed]['value'])
    log(f'check finch/goa/{seed}: AUPRC {value:.10f}, scorer {want:.10f}')
    if abs(value - want) > 1e-6:
        raise ValueError('the official labels do not reproduce the scorer')
    reps = {ts: np.mean([replicate_auprc(scores_of(p), y_official, keep_all, cluster_of, draws) for p in files['finch', ts].values()], 0)
            for ts in ('archipelago', 'hq_random')}
    low, high = np.percentile(reps['archipelago'] - reps['hq_random'], [2.5, 97.5])
    row = next(r for r in csv.DictReader(open(results / 'main' / 'contrasts.tsv', newline=''), delimiter='\t')
               if r['metric'] == 'auprc' and (r['method'], r['training_set'], r['minus_method'], r['minus_training_set'])
               == ('finch', 'archipelago', 'finch', 'hq_random'))
    log(f"check interval finch archipelago - hq_random: [{low:.6f}, {high:.6f}], scorer [{float(row['ci_low']):.6f}, "
        f"{float(row['ci_high']):.6f}]")
    if max(abs(low - float(row['ci_low'])), abs(high - float(row['ci_high']))) > 1e-5:
        raise ValueError('the bootstrap does not reproduce the scorer interval')

    # ---- the versions: (labels, pairs kept)
    def rows_of(proteins):
        rows = torch.zeros(len(ids), dtype=torch.bool, device=dev)
        rows[[row_of[h] for h in proteins]] = True
        return rows[:, None].expand_as(keep_all).contiguous()

    y_exp = tensor(closure(pairs['exp']))
    subsets = released_subsets(root)
    versions = {'experimental': (y_exp, y_exp.any(1)[:, None].expand_as(keep_all).contiguous()),
                'direct_negatives': (tensor(direct), keep_all), 'direct_masked': (tensor(direct), tensor(~inferred)),
                'masked': (tensor(base), tensor(~masked)), 'reviewed': (y_official, keep_all),
                'all_goa': (tensor(closure(pairs['all_goa'])), keep_all),
                'no_hit_0.001': (y_official, rows_of(subsets[1e-3])), 'no_hit_10': (y_official, rows_of(subsets[10.0]))}

    # ---- 1. the direct-evidence versions and the computational annotations as negatives: seed means and paired intervals
    tables = {**versions, 'masked_negatives': (tensor(base), keep_all)}  # (not a column of the grid)
    for name in TABLES:
        y, keep = tables[name]
        point, replicates, start = {}, {}, time.time()
        outs[name].mkdir(parents=True)
        with open(outs[name] / 'models.tsv', 'w', newline='\n') as handle:
            handle.write('method\ttraining_set\tseed\tmetric\tvalue\n')
            for key in BOOTSTRAP:
                values, per_seed_reps = [], []
                for seed, path in files[key].items():
                    s = scores_of(path)
                    values.append(auprc(s, y, keep))
                    per_seed_reps.append(replicate_auprc(s, y, keep, cluster_of, draws))
                    handle.write(f'{key[0]}\t{key[1]}\t{seed}\tauprc\t{values[-1]}\n')
                point[key], replicates[key] = float(np.mean(values)), np.mean(per_seed_reps, 0)
                handle.write(f'{key[0]}\t{key[1]}\tmean\tauprc\t{point[key]}\n')
                log(f'{time.time() - start:6.0f} s', name, *key, round(point[key], 4))
        with open(outs[name] / 'contrasts.tsv', 'w', newline='\n') as handle:
            handle.write('metric\tmethod\ttraining_set\tminus_method\tminus_training_set\tdifference\tci_low\tci_high\n')
            for x in BOOTSTRAP:
                for z in BOOTSTRAP:
                    if x != z and (x[0] == z[0] or x[1] == z[1]):
                        low, high = np.percentile(replicates[x] - replicates[z], [2.5, 97.5])
                        handle.write(f'auprc\t{x[0]}\t{x[1]}\t{z[0]}\t{z[1]}\t{point[x] - point[z]}\t{low}\t{high}\n')

    # ---- 2. the grid, twice: results/label_grid with the two no-hit columns on the reviewed labels, and the paper's grid,
    # results_masked/label_grid, with them on the main test labels (every other column is the same in both)
    main_versions = {**versions, 'no_hit_0.001': (tensor(base), tensor(~masked) & rows_of(subsets[1e-3])),
                     'no_hit_10': (tensor(base), tensor(~masked) & rows_of(subsets[10.0]))}
    grids = {outs['grid']: (versions, [(v, results / f) for v, f in REFERENCES]),
             outs['main_grid']: (main_versions, [(v, results / f) for v, f in REFERENCES if not v.startswith('no_hit')]
                                 + [(v, results_masked / 'no_hit' / f'E{c}') for v, c in (('no_hit_0.001', '0.001'),
                                                                                          ('no_hit_10', '10'))])}
    point, start = {out: {} for out in grids}, time.time()
    for key in GRID_ROWS:
        per_seed = collections.defaultdict(list)
        for path in files[key].values():
            s = scores_of(path)
            for name in GRID_VERSIONS:
                per_seed[name].append(auprc(s, *versions[name]))
                if main_versions[name] is not versions[name]:
                    per_seed['main', name].append(auprc(s, *main_versions[name]))
        for name in GRID_VERSIONS:
            point[outs['grid']][(*key, name)] = float(np.mean(per_seed[name]))
            point[outs['main_grid']][(*key, name)] = float(np.mean(per_seed[('main', name) if ('main', name) in per_seed
                                                                            else name]))
        log(f'{time.time() - start:6.0f} s', *key, len(files[key]), 'runs',
            {n: round(point[outs['main_grid']][(*key, n)], 4) for n in GRID_VERSIONS}, '(main labels in the no-hit columns)')

    for out, (labels, references) in grids.items():
        out.mkdir(parents=True)
        with open(out / 'versions.tsv', 'w', newline='\n') as handle:
            handle.write('version\tproteins\tclusters\tpositive_pairs\tleft_out_pairs\tterms_with_positive\n')
            for name, (y, keep) in labels.items():
                rows = keep.any(1)
                handle.write(f'{name}\t{int(rows.sum())}\t{len(torch.unique(cluster_of[rows]))}\t{int((y & keep).sum())}\t'
                             f'{int((~keep & rows[:, None]).sum())}\t{int((y & keep).any(0).sum())}\n')
        log((out / 'versions.tsv').read_text())
        with open(out / 'grid.tsv', 'w', newline='\n') as handle:
            handle.write('method\ttraining_set\truns\t' + '\t'.join(GRID_VERSIONS) + '\n')
            for key in GRID_ROWS:
                handle.write(f'{key[0]}\t{key[1]}\t{len(files[key])}\t'
                             + '\t'.join(repr(point[out][(*key, n)]) for n in GRID_VERSIONS) + '\n')

        checked, worst = [], 0.0
        for version, folder in references:
            path = folder / ('models.tsv' if version in TABLES else 'metrics.tsv')
            if not path.exists():
                raise FileNotFoundError(f'{path} is missing: the grid is checked against it (run that scoring stage first)')
            for r in csv.DictReader(open(path, newline=''), delimiter='\t'):
                key = (r['method'], r['training_set'], version)
                if r['seed'] == 'mean' and r['metric'] == 'auprc' and r.get('split', 'test') == 'test' and key in point[out]:
                    checked.append((*key, point[out][key], float(r['value']), folder.relative_to(root).as_posix()))
                    worst = max(worst, abs(point[out][key] - float(r['value'])))
        with open(out / 'check.tsv', 'w', newline='\n') as handle:
            handle.write('method\ttraining_set\tversion\tgrid\treference\tdifference\treference_folder\n')
            handle.writelines(f'{m}\t{t}\t{v}\t{x!r}\t{z!r}\t{x - z!r}\t{f}\n' for m, t, v, x, z, f in checked)
        log(f'{out}: checked {len(checked)} grid values against the scored results '
            f'{dict(collections.Counter(v for _, _, v, *_ in checked))}; largest difference {worst:.2e}')
        if worst > 1e-6:
            raise ValueError(f'a grid value differs from its reference ({out / "check.tsv"})')
    log('wrote', *outs.values())

if __name__ == '__main__':
    main()
