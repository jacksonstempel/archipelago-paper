"""The data comparisons on the terms each original training set covers (Supplementary Note 10): the paper's scorer
(archipelago.scorer), run unchanged on the fits of the main run, with every metric restricted to the GO terms that have at
least one positive in one training set.

    python -m retrain.scoring.shared_terms [--data DIR] [--labels reviewed|masked] [--device cuda|cpu]

1. For `all` (every term: the scorer as it is) and for each published method's original data (original_deepgose,
   original_transfew, original_dpfunc), results/shared_terms/<set>/: the scorer's results folder, its contrasts.tsv
   (retrain.scoring.score.contrasts) and terms.json (the covered terms and, per split, the label pairs and proteins on the
   terms the set lacks). The restriction is made inside the scorer's Benchmark, for the test and the development proteins
   alike (so F1's development threshold is chosen on the same terms): the scores are propagated over the whole vocabulary
   as usual, then the labels, the information accretion and the propagated scores keep only the covered terms' columns. A
   protein left without a positive is not evaluated (the scorer needs a label per protein); the bootstrap keeps the
   scorer's draws of the test clusters, so every replicate is the paper's (a cluster left without proteins adds nothing to
   it).
2. The check: every row of results/shared_terms/all/metrics.tsv (test and development, every seed, mean and ensemble;
   value, threshold, interval, seed SD) and every replicate array equals results/main to 1e-9.
3. results/shared_terms/summary.tsv and summary.md: for each original set, the paper's data comparisons over every term
   and over the terms the set covers (seed means, paired 95% cluster-bootstrap intervals): the published method from its
   original data to HQ to ArchipelaGO, Finch from the set to ArchipelaGO, and the spread of the three published methods
   on ArchipelaGO (best minus worst; for Smin, lower is better).

--labels masked: the main test labels (the pairs the reviewed GOA assertions add are neither positive nor negative; those on the covered\nterms stay out), from results_masked/main into results_masked/shared_terms.\n\nNeeds results/main (retrain.scoring.score main). About a minute per set on one GPU.
"""
import argparse
import csv
import json
import os
import subprocess
import sys
import traceback
from pathlib import Path

import numpy as np

from retrain.common import add_argument, data_root, load, log
from retrain.scoring.score import RESULTS, contrasts, dataset, fits, results_dir, scorer_labels, write_config

SETS = {'deepgose': 'original_deepgose', 'transfew': 'original_transfew', 'dpfunc': 'original_dpfunc'}
NAMES = {'deepgose': 'DeepGO-SE', 'transfew': 'TransFew', 'dpfunc': 'DPFunc', 'finch': 'Finch'}
METRICS = ('auprc', 'fmax', 'smin')


def restricted(benchmark, columns, stats):
    """The scorer's Benchmark class over `columns` and the proteins with a positive among them; each split's counts go into
    `stats`."""
    import torch

    class Restricted(benchmark):
        def __init__(self, labels, ia=None, clusters=None, draws=None, device='cpu', mask=None, select=None):
            y = torch.as_tensor(labels, device=device).bool()
            self._columns = torch.as_tensor(columns, device=device)
            kept = y.index_select(1, self._columns)
            off = y.clone()
            off[:, self._columns] = False
            rows = torch.nonzero(kept.any(1)).squeeze(1)
            split = 'development' if clusters is None else 'test'
            stats[split] = dict(proteins=len(y), positive_pairs=int(y.sum()), positive_pairs_on_uncovered_terms=int(off.sum()),
                                share_on_uncovered_terms=float(off.sum()) / float(y.sum()),
                                uncovered_terms_with_a_positive=int(off.any(0).sum()),
                                proteins_with_a_positive_on_uncovered_terms=int(off.any(1).sum()),
                                proteins_evaluated=len(rows), positive_pairs_evaluated=int(kept.sum()))
            if mask is not None:
                mask = torch.as_tensor(mask, device=device).bool().index_select(1, self._columns)
            if len(rows) < len(y):
                kept = kept.index_select(0, rows)
                mask = None if mask is None else mask.index_select(0, rows)
                at = rows.cpu().numpy()
                select = at if select is None else np.asarray(select)[at]
                if clusters is not None:
                    clusters = np.asarray(clusters)[at]
            if clusters is not None:
                stats[split].update(clusters=int(np.asarray(draws).shape[1]), clusters_with_proteins=len(np.unique(clusters)))
                # clusters after the last one with proteins: no pairs in any replicate, so their columns of the draws go
                draws = np.asarray(draws)[:, :int(np.max(clusters)) + 1]
            log(split, stats[split])
            ia = None if ia is None else np.asarray(ia.cpu() if torch.is_tensor(ia) else ia)[columns]
            super().__init__(kept, ia, clusters, draws, device, mask=mask, select=select)

        def statistics(self, scores):  # scores: propagated over the whole vocabulary
            return super().statistics(torch.as_tensor(scores, device=self.device).index_select(1, self._columns))

    return Restricted


def score_set(root, name, config, device, labels='reviewed'):
    """In this process: the scorer on `config`, restricted to the terms of `name` (or every term), then terms.json."""
    from archipelago import dataset as released, evaluate as ev, scorer
    vocabulary = json.loads((root / 'datasets' / 'ontology' / 'vocabulary.json').read_text(encoding='utf-8'))
    data = config.parent / json.loads(config.read_text(encoding='utf-8'))['dataset']
    if released.terms(data)[0] != vocabulary:
        raise ValueError('the scorer vocabulary differs from the ontology dataset')
    if name == 'all':
        covered, source = np.ones(len(vocabulary), dtype=bool), None
    else:
        d = load(root, name, sequences=False)
        covered = np.asarray((d.labels > 0).sum(0)).ravel() > 0
        source = dict(dataset=name, proteins=len(d.seq_ids))
    columns = np.flatnonzero(covered)
    log(f'{name}: {len(columns)} of {len(vocabulary)} terms')
    stats = {}
    ev.Benchmark = restricted(ev.Benchmark, columns, stats)  # the scorer's main() takes Benchmark from archipelago.evaluate
    out = results_dir(root, labels) / 'shared_terms' / name
    scorer.main([str(config), '--out', str(out), '--cache', str(root / 'cache' / 'scoring'), '--labels', scorer_labels(labels)]
                + (['--device', device] if device else []))
    record = dict(source=source, terms=len(columns), vocabulary=len(vocabulary), splits=stats,
                  covered_terms=[vocabulary[i] for i in columns], uncovered_terms=[t for t, c in zip(vocabulary, covered) if not c])
    (out / 'terms.json').write_text(json.dumps(record, indent=1) + '\n', encoding='utf-8')


def read(path):
    with open(path, encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def check(results):
    """results/shared_terms/all against results/main: every metrics row and replicate array to 1e-9."""
    w, db = results / 'shared_terms', results / 'main'
    key = ('method', 'training_set', 'seed', 'split', 'metric')
    mine = {tuple(r[k] for k in key): r for r in read(w / 'all' / 'metrics.tsv')}
    official = {tuple(r[k] for k in key): r for r in read(db / 'metrics.tsv')}
    if set(mine) != set(official):
        raise ValueError(f'rows differ: {sorted(set(mine) ^ set(official))[:5]}')
    worst = 0.0
    for k, r in official.items():
        for c in ('value', 'threshold', 'ci_low', 'ci_high', 'seed_sd'):
            if (r[c] in ('', 'nan')) != (mine[k][c] in ('', 'nan')):
                raise ValueError(f'{k} {c}: {r[c]!r} vs {mine[k][c]!r}')
            if r[c] not in ('', 'nan'):
                worst = max(worst, abs(float(r[c]) - float(mine[k][c])))
    za, zb = np.load(db / 'replicates.npz'), np.load(w / 'all' / 'replicates.npz')
    if set(za.files) != set(zb.files):
        raise ValueError('replicate arrays differ')
    worst_rep = max(float(np.abs(za[k] - zb[k]).max()) for k in za.files if k != 'metrics')
    text = (f'shared_terms/all against main: {len(official)} rows (test and development, every seed, mean and ensemble), '
            f'largest difference {worst:.1e}; {len(za.files) - 1} replicate arrays, largest difference {worst_rep:.1e}')
    log(text)
    if worst > 1e-9 or worst_rep > 1e-9:
        raise ValueError('the full-vocabulary scoring does not reproduce results/main')
    return text


def summary(w, check_text):
    """w/summary.tsv and w/summary.md from the folders w/all and w/<original set>."""
    point, ci, terms = {}, {}, {}
    for name in ['all', *SETS.values()]:
        terms[name] = json.loads((w / name / 'terms.json').read_text(encoding='utf-8'))
        for r in read(w / name / 'metrics.tsv'):
            if r['seed'] == 'mean' and r['split'] == 'test':
                point[name, r['method'], r['training_set'], r['metric']] = float(r['value'])
        for r in read(w / name / 'contrasts.tsv'):
            ci[name, r['metric'], r['method'], r['training_set'], r['minus_method'], r['minus_training_set']] = \
                (float(r['difference']), float(r['ci_low']), float(r['ci_high']))
    rows = []

    def add(comparison, name, a, b):
        for metric in METRICS:
            d, lo, hi = ci[(name, metric, *a, *b)]
            va, vb = point[(name, *a, metric)], point[(name, *b, metric)]
            if abs(d - (va - vb)) > 1e-12:
                raise ValueError(f'{comparison} {name} {metric}: contrast {d} vs points {va - vb}')
            rows.append(dict(comparison=comparison, terms=name, n_terms=terms[name]['terms'], metric=metric,
                             method=a[0], training_set=a[1], minus_method=b[0], minus_training_set=b[1], value=va,
                             minus_value=vb, difference=d, ci_low=lo, ci_high=hi))

    for method, own in SETS.items():
        for name in ('all', own):
            add(f'{NAMES[method]}, original to HQ', name, (method, 'hq'), (method, 'original'))
            add(f'{NAMES[method]}, HQ to ArchipelaGO', name, (method, 'archipelago'), (method, 'hq'))
            add(f'{NAMES[method]}, original to ArchipelaGO', name, (method, 'archipelago'), (method, 'original'))
            add(f"Finch, {NAMES[method]}'s original set to ArchipelaGO", name, ('finch', 'archipelago'), ('finch', own))
    for name in ('all', *SETS.values()):
        for metric in METRICS:  # spread: best minus worst published method on ArchipelaGO, by the point estimates
            value = {m: point[name, m, 'archipelago', metric] for m in SETS}
            best, worst = sorted(value, key=value.get, reverse=metric != 'smin')[::2]
            d, lo, hi = ci[name, metric, best, 'archipelago', worst, 'archipelago']
            rows.append(dict(comparison='Spread of the published methods on ArchipelaGO', terms=name, n_terms=terms[name]['terms'],
                             metric=metric, method=best, training_set='archipelago', minus_method=worst,
                             minus_training_set='archipelago', value=value[best], minus_value=value[worst], difference=d,
                             ci_low=lo, ci_high=hi))
    with open(w / 'summary.tsv', 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

    fmt = {'auprc': '{:+.3f}', 'fmax': '{:+.3f}', 'smin': '{:+.2f}'}
    val = {'auprc': '{:.3f}', 'fmax': '{:.3f}', 'smin': '{:.2f}'}
    label = lambda name: 'every term' if name == 'all' else f'shared terms ({name})'  # noqa: E731
    pick = lambda c, n, m: next(x for x in rows if x['comparison'] == c and x['terms'] == n and x['metric'] == m)  # noqa: E731
    shown = {'auprc': 'AUPRC', 'fmax': 'Fmax', 'smin': 'Smin'}
    lines = ['# Data gains on the terms each original training set covers', '', check_text, '',
             '## Test label pairs on terms each original set lacks', '',
             '| Original set | Terms covered | Test pairs on uncovered terms | Share | Test proteins with such a pair | '
             'Uncovered terms with a test positive | Test proteins left without a positive | Development pairs on uncovered terms |',
             '|---|---:|---:|---:|---:|---:|---:|---:|']
    for method, own in SETS.items():
        t, s, dv = terms[own], terms[own]['splits']['test'], terms[own]['splits']['development']
        lines.append(f"| {NAMES[method]}'s ({own}) | {t['terms']:,} | {s['positive_pairs_on_uncovered_terms']:,} of "
                     f"{s['positive_pairs']:,} | {100 * s['share_on_uncovered_terms']:.2f}% | "
                     f"{s['proteins_with_a_positive_on_uncovered_terms']:,} | {s['uncovered_terms_with_a_positive']:,} | "
                     f"{s['proteins'] - s['proteins_evaluated']} | {dv['positive_pairs_on_uncovered_terms']:,} of "
                     f"{dv['positive_pairs']:,} ({100 * dv['share_on_uncovered_terms']:.2f}%) |")
    lines += ['', '## Comparisons: every term and shared terms (test; seed means; paired 95% cluster-bootstrap intervals)', '',
              'Smin: lower is better, so a gain is negative.', '',
              '| Comparison | Terms | AUPRC from → to | AUPRC difference | Fmax difference | Smin difference |',
              '|---|---|---|---|---|---|']
    for r in rows:
        if r['metric'] != 'auprc':
            continue
        cells, which = [], []
        for metric in METRICS:
            q = pick(r['comparison'], r['terms'], metric)
            cells.append(f"{fmt[metric].format(q['difference'])} [{fmt[metric].format(q['ci_low'])}, "
                         f"{fmt[metric].format(q['ci_high'])}]")
            which.append(f"{shown[metric]} {NAMES[q['method']]} − {NAMES[q['minus_method']]}")
        what = r['comparison'] + (f" (best − worst: {'; '.join(which)})" if 'Spread' in r['comparison'] else '')
        lines.append(f"| {what} | {label(r['terms'])}, {r['n_terms']:,} | {val['auprc'].format(r['minus_value'])} → "
                     f"{val['auprc'].format(r['value'])} | " + ' | '.join(cells) + ' |')
    lines += ['', '## Data gain (original to ArchipelaGO) against the spread on ArchipelaGO, on the same terms', '',
              '| Method | Terms | Gain AUPRC | Spread AUPRC | Ratio | Gain Fmax | Spread Fmax | Ratio | Gain Smin | Spread Smin | Ratio |',
              '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for method, own in SETS.items():
        for name in ('all', own):
            g = {m: pick(f'{NAMES[method]}, original to ArchipelaGO', name, m)['difference'] for m in METRICS}
            s = {m: pick('Spread of the published methods on ArchipelaGO', name, m)['difference'] for m in METRICS}
            lines.append(f"| {NAMES[method]} | {label(name)} | " + ' | '.join(
                f"{fmt[m].format(g[m])} | {fmt[m].format(s[m])} | {g[m] / s[m]:.1f}" for m in METRICS) + ' |')
    lines += ['', '## Test values per model (seed means)', '', '| Model, training set | Terms | AUPRC | Fmax | Smin |',
              '|---|---|---:|---:|---:|']
    models = [(m, s) for m in SETS for s in ('original', 'hq', 'archipelago')] + \
             [('finch', s) for s in (*SETS.values(), 'hq', 'archipelago')]
    for m, s in models:
        for name in ['all', *SETS.values()]:
            lines.append(f"| {NAMES[m]}, {s} | {label(name)} | " +
                         ' | '.join(val[k].format(point[name, m, s, k]) for k in METRICS) + ' |')
    (w / 'summary.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    log(f'wrote {w / "summary.tsv"} and summary.md')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--device', choices=['cuda', 'cpu'], help="the scorer's device (default: cuda if available)")
    parser.add_argument('--set', help=argparse.SUPPRESS)  # internal: score one set in this process
    parser.add_argument('--labels', choices=list(RESULTS), default='reviewed', help='the version of the main test labels')
    add_argument(parser)
    a = parser.parse_args(argv)
    root = data_root(a.data)
    # masked: the predictions of the reviewed scoring of the main run
    config = root / 'scoring' / ('shared_terms.json' if a.labels == 'reviewed' else 'main.json')
    if a.set:  # as the scorer's own entry point: exit at once, without waiting for its reader threads
        try:
            score_set(root, a.set, config, a.device, a.labels)
        except BaseException:
            traceback.print_exc()
            sys.stderr.flush()
            os._exit(1)
        sys.stdout.flush()
        os._exit(0)
    results = results_dir(root, a.labels)
    if not (results / 'main' / 'replicates.npz').is_file():
        raise FileNotFoundError(f'{results / "main"}: run retrain.scoring.score main (with --labels {a.labels}) first')
    (results / 'shared_terms').mkdir(parents=True, exist_ok=False)
    if a.labels == 'reviewed':
        write_config(config, dataset(root), fits(root, 'main'))
    for name in ('all', *SETS.values()):
        command = [sys.executable, '-m', 'retrain.scoring.shared_terms', '--set', name, '--labels', a.labels, '--data', str(root)]
        log('RUN', ' '.join(command + (['--device', a.device] if a.device else [])))
        subprocess.run(command + (['--device', a.device] if a.device else []), check=True)
        contrasts(results / 'shared_terms' / name)
    summary(results / 'shared_terms', check(results))


if __name__ == '__main__':
    main()
