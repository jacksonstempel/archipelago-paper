"""Score a set of fits with the paper's scorer (archipelago.scorer) and write their paired contrasts.

    python -m retrain.scoring.score RUN [--data DIR] [--device cuda|cpu]

RUN is one of
  main                every model on every training set of the main comparisons (Table 1, Figs. 2 and 4b, SI tables)
  selection           Finch on HQ plus 150,000 proteins chosen by each rule of Fig. 3b, and with 300,000 / 600,000 picks
  ladder              the models of the Finch recipe ablation (Fig. 4a; configs/finch.json, 'ladder'), Finch's variants
                      ('variants'), Finch and the published methods, all trained on ArchipelaGO
  explabels_training  every model trained on the HQ proteins with their experimental labels only (hq_explabels), beside
                      its original data, HQ and ArchipelaGO (Supplementary Information)

Steps:
1. The scorer's data bundle, <root>/scoring/data/ (the format archipelago.scorer reads): vocabulary.txt, go-basic.obo,
   benchmark/{development,test}.tsv and their _labels.tsv, and archipelago/{proteins,labels}.tsv holding the HQ proteins
   (the scorer estimates the information accretion of Smin from them). Written from datasets/ (ontology, development,
   test, hq).
2. The scorer's config, <root>/scoring/RUN.json: every fit runs/<method>/<training set>/seed<k>/ of the run (every seed
   folder found; a fit with record.json and both score files). TransFew's further checkpoints (side files bcebest,
   auprcbest) are scored as the methods transfew_bcebest and transfew_auprcbest. A missing fit stops the export.
3. archipelago-score CONFIG --out <root>/results/RUN (which must not exist yet).
4. <root>/results/RUN/contrasts.tsv: for every metric (auprc, fmax, f1, smin), every pair of training sets of one model and
   every pair of models on one training set: the difference of the seed means and its paired 95% cluster-bootstrap
   interval (2.5th and 97.5th percentiles of the replicate differences in replicates.npz).
"""
import argparse
import csv
import itertools
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

from archipelago.benchmark import write_bundle
from retrain.common import add_argument, data_root, load, log, score_file
from retrain.finch import config as finch_config

PUBLISHED = ('deepgose', 'dpfunc', 'transfew')
BASELINES = ('naive', 'diamond_knn', 'embedding_knn')
SIDE = {'transfew': ('bcebest', 'auprcbest')}  # further checkpoints, scored as <method>_<name>
RUNS = {  # run: {method: training sets}
    'main': {'finch': ('goa', 'hq', 'hq_random', 'archipelago'),
             **{m: ('original', 'goa', 'hq', 'archipelago') for m in PUBLISHED},
             **{m: ('goa', 'hq', 'archipelago') for m in BASELINES}},
    'selection': {'finch': ('hq', 'hq_random', 'pool_random', 'filter_random', 'coverage_unfiltered', 'archipelago',
                            'archipelago_300k', 'archipelago_600k')},
    'ladder': {m: ('archipelago',) for m in (*finch_config()['ladder'], *finch_config()['variants'], 'finch', *PUBLISHED)},
    'explabels_training': {'finch': ('hq_explabels', 'hq', 'archipelago'),
                           **{m: ('original', 'hq_explabels', 'hq', 'archipelago') for m in PUBLISHED},
                           **{m: ('hq_explabels', 'hq', 'archipelago') for m in BASELINES}},
}
CONTRAST_METRICS = ('auprc', 'fmax', 'f1', 'smin')


def bundle(root):
    """Write <root>/scoring/data/ from the dataset folders (archipelago.benchmark.write_bundle); returns its path."""
    data = Path(root) / 'scoring' / 'data'
    shutil.rmtree(data, ignore_errors=True)
    ontology = Path(root) / 'datasets' / 'ontology'
    vocabulary = json.loads((ontology / 'vocabulary.json').read_text(encoding='utf-8'))
    sets = {}
    for name in ('development', 'test', 'hq'):
        ds = load(root, name, scoring=True, sequences=False)
        rows = [dict(protein_id=h, accessions=a, **({'cluster': ds.columns['cluster'][i]} if name != 'hq' else {}))
                for i, (h, a) in enumerate(zip(ds.seq_ids, ds.accessions))]
        y = ds.labels.tocoo()
        sets[name] = (rows, [(ds.seq_ids[r], vocabulary[c]) for r, c in zip(y.row.tolist(), y.col.tolist())])
    return write_bundle(data, vocabulary, ontology / 'go-basic.obo', sets)


def fits(root, run):
    """[(method as scored, training set, {seed: {'development': path, 'test': path}})] of a run; fails on a missing fit."""
    out = []
    for method, training_sets in RUNS[run].items():
        for side in (None, *SIDE.get(method, ())):
            for training_set in training_sets:
                folder = Path(root) / 'runs' / method / training_set
                seeds = {}
                for seed_dir in sorted(folder.glob('seed*'), key=lambda p: int(p.name[4:])):
                    files = {split: score_file(seed_dir, prefix, side) for split, prefix in
                             (('development', 'dev'), ('test', 'test'))}
                    if not (seed_dir / 'record.json').is_file() or not all(p.is_file() for p in files.values()):
                        raise FileNotFoundError(f'{seed_dir}: incomplete fit (record.json or a score file missing)')
                    seeds[seed_dir.name[4:]] = files
                if not seeds:
                    raise FileNotFoundError(f'{folder}: no fit of {method} on {training_set}')
                out.append((method + (f'_{side}' if side else ''), training_set, seeds))
    return out


def write_config(path, data, models, test_files=None):
    """A scorer config at `path`: paths relative to it. test_files optionally replaces each test file
    ({(method, training set, seed): path})."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rel = lambda p: os.path.relpath(p, path.parent).replace(os.sep, '/')  # noqa: E731
    listed = [dict(method=m, training_set=t, seeds={
        seed: {'development': rel(files['development']),
               'test': rel((test_files or {}).get((m, t, seed), files['test']))} for seed, files in seeds.items()})
        for m, t, seeds in models]
    path.write_text(json.dumps(dict(data=rel(data), models=listed), indent=1) + '\n')
    return path


def test_subset(root, data, folder, keep, test_labels, models):
    """A scorer bundle and config in `folder` that score `models` on the test proteins `keep` (a set of seq_ids) with the
    test label rows `test_labels` ([(seq_id, go_id)]); everything else as in `data`. Each test prediction is a symbolic
    link with a .ids.txt row list (the full test order), so the scorer takes the subset's rows by identifier."""
    folder = Path(folder)
    if folder.exists():
        raise FileExistsError(f'{folder} exists')
    subset = folder / 'data'
    shutil.copytree(data, subset)
    with open(data / 'benchmark' / 'test.tsv', newline='') as handle:
        table = list(csv.reader(handle, delimiter='\t'))
    with open(subset / 'benchmark' / 'test.tsv', 'w', newline='\n') as handle:
        handle.writelines('\t'.join(row) + '\n' for row in [table[0]] + [r for r in table[1:] if r[0] in keep])
    with open(subset / 'benchmark' / 'test_labels.tsv', 'w', newline='\n') as handle:
        handle.write('protein_id\tgo_id\n')
        handle.writelines(f'{h}\t{g}\n' for h, g in test_labels if h in keep)
    order = '\n'.join(r[0] for r in table[1:]) + '\n'
    links = {}
    for method, training_set, seeds in models:
        for seed, files in seeds.items():
            link = folder / 'predictions' / method / training_set / f'seed{seed}' / 'test_scores.npy'
            link.parent.mkdir(parents=True)
            os.symlink(Path(files['test']).resolve(), link)
            with open(link.with_suffix('.ids.txt'), 'w', newline='\n') as handle:
                handle.write(order)
            links[method, training_set, seed] = link
    return write_config(folder / 'config.json', subset, models, links)


def run_scorer(root, config, out, device=None):
    """The scorer (archipelago.scorer) on `config`, results in `out` (which must not exist)."""
    command = [sys.executable, '-m', 'archipelago.scorer', str(config), '--out', str(out),
               '--cache', str(Path(root) / 'cache' / 'scoring')]
    if device:
        command += ['--device', device]
    log('RUN', ' '.join(command))
    subprocess.run(command, check=True)


def contrasts(results):
    """results/contrasts.tsv: paired differences of seed means (see the module docstring); returns the rows."""
    from archipelago.evaluate import paired_difference
    results = Path(results)
    replicates = np.load(results / 'replicates.npz')
    names = [str(m) for m in replicates['metrics']]
    with open(results / 'metrics.tsv', newline='') as handle:
        point = {(r['method'], r['training_set'], r['metric']): float(r['value'])
                 for r in csv.DictReader(handle, delimiter='\t') if r['seed'] == 'mean' and r['split'] == 'test'}
    models = [tuple(k.split('/')[:2]) for k in replicates.files if k.endswith('/mean')]
    pairs = [((m, a), (m, b)) for m in sorted({m for m, _ in models})
             for a, b in itertools.permutations([s for mm, s in models if mm == m], 2)]
    pairs += [((m1, s), (m2, s)) for s in sorted({s for _, s in models})
              for m1, m2 in itertools.permutations([m for m, ss in models if ss == s], 2)]
    rows = []
    for metric in CONTRAST_METRICS:
        column = names.index(metric)
        for (m1, s1), (m2, s2) in pairs:
            difference, (low, high) = paired_difference(
                (point[m1, s1, metric], replicates[f'{m1}/{s1}/mean'][:, column]),
                (point[m2, s2, metric], replicates[f'{m2}/{s2}/mean'][:, column]))
            rows.append(dict(metric=metric, method=m1, training_set=s1, minus_method=m2, minus_training_set=s2,
                             difference=difference, ci_low=float(low), ci_high=float(high)))
    with open(results / 'contrasts.tsv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    log(f'wrote {results / "contrasts.tsv"}: {len(rows):,} contrasts')
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('run', choices=list(RUNS))
    parser.add_argument('--device', choices=['cuda', 'cpu'], help="the scorer's device (default: cuda if available)")
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    data = bundle(root)
    config = write_config(root / 'scoring' / f'{a.run}.json', data, fits(root, a.run))
    out = root / 'results' / a.run
    run_scorer(root, config, out, a.device)
    contrasts(out)


if __name__ == '__main__':
    main()
