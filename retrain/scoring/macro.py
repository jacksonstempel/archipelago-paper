"""Term-centric (macro) AUPRC: the average precision of each GO term over the test proteins, averaged over the terms with at
least one test positive (Supplementary Note 11 and Table 15; the paper's scorer has no macro metric).

    python -m retrain.scoring.macro [--data DIR] [--labels reviewed|masked] [--device cuda|cpu]

Every model and seed of the scorer configs of the main comparison and of the training run with experimental labels only
(<root>/scoring/{main,explabels_training}.json, written by retrain.scoring.score; a model listed in both, with the same test
files, once): the raw test scores (rows as the scorer reads them) propagated with the scorer's Propagation, then
- micro-AUPRC and its cluster-bootstrap replicates from the scorer's own Benchmark. Check: every seed's value and
  replicates equal its results folder's (results/main, results/explabels_training) to 1e-9, or stop.
- macro-AUPRC: for each term with a test positive, the exact average precision over the test proteins (tied scores grouped,
  as the scorer's micro-AUPRC), averaged over those terms.
  Intervals (paired; 10,000 replicates): a bootstrap over terms, the terms drawn with replacement (seed 20261015, the
  scorer's), each term's average precision fixed. The paper's cluster bootstrap is computed too, as a diagnostic: in a
  replicate a protein counts as many times as its cluster was drawn and the average is over the terms with a drawn positive.
  It is not used, because about a fifth of the terms have no drawn positive in a replicate, so each replicate averages a
  different set of terms, without the rarest; its percentile intervals then miss the point estimates (models.tsv,
  macro_cluster_* columns; contrasts.tsv, bootstrap `clusters (diagnostic, not used)`).
  Checks on the first run: the point value equals the mean of the scorer's micro_auprc applied to each term's column, and
  two cluster replicates equal that reference on the data set with the proteins repeated.
--labels masked: the masked test labels (every metric leaves the reviewed GOA assertions' pairs out: neither positive nor negative, in the
micro-AUPRC and in each term's average precision); the databanks are results_masked/{main,explabels_training}, the output
results_masked/macro.
Seed means; a model's replicates are the means of its seeds' (as the scorer). Writes results/macro/models.tsv, seeds.tsv,
contrasts.tsv (the paper's comparisons; figures/si.py reads it), contrasts_all.tsv (every pair of sets within a model and
of models within a set; micro over clusters, macro over terms), replicates.npz, per_term.npz (each model's per-term
average precision, seed mean) and check.txt. Needs results/main and results/explabels_training (retrain.scoring.score);
about a minute on one GPU.
"""
import argparse
import csv
import hashlib
import io
import itertools
import json
import os
import time
from pathlib import Path

import numpy as np
import torch

from archipelago import dataset, evaluate as ev, scorer
from archipelago.go import Propagation
from retrain.paths import add_argument, data_root
from retrain.scoring.score import RESULTS, results_dir, scorer_labels


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--device', default='cuda', help='cuda (the paper) or cpu')
    parser.add_argument('--labels', choices=list(RESULTS), default='reviewed', help='the version of the main test labels')
    add_argument(parser)
    a = parser.parse_args(argv)
    root, device, labels = data_root(a.data), a.device, a.labels
    out = results_dir(root, labels) / 'macro'
    pairs = [(root / 'scoring' / f'{run}.json', results_dir(root, labels) / run)
             for run in ('main', 'explabels_training')]
    CACHE = root / 'cache' / 'scoring'
    dev = torch.device(device)
    torch.backends.cuda.matmul.allow_tf32 = False  # the count products below must be exact (integers below 2**24 in float32)
    torch.set_float32_matmul_precision('highest')
    out.mkdir(parents=True, exist_ok=False)
    log = open(out / 'check.txt', 'w')

    def say(*args):
        text = ' '.join(str(a) for a in args)
        print(text, flush=True)
        log.write(text + '\n')
        log.flush()

    # ---- models of every config (a model in several configs once, if its test files are the same)
    models, data = {}, set()
    for config, databank in pairs:
        c = json.load(open(config))
        base = Path(config).parent
        data.add(os.path.realpath(base / c['dataset']))
        for m in c['models']:
            key = (m['method'], m['training_set'])
            files = {seed: os.path.realpath(base / f['test']) for seed, f in m['seeds'].items()}
            if key in models:
                if models[key]['files'] != files:
                    raise ValueError(f'{key} has different test files in {config}')
                continue
            models[key] = dict(files=files, databank=Path(databank))
    if len(data) != 1:
        raise ValueError(f'the configs use different dataset folders: {data}')
    data = Path(data.pop())
    say(f'{len(models)} models, {sum(len(m["files"]) for m in models.values())} runs; data {data}')

    # ---- the scorer's inputs (its cache): test labels, clusters and draws; information accretion and parents from go_terms
    scorer.ev = ev  # prepared_inputs uses the module the scorer's main() imports
    vocabulary, parents, ia = dataset.terms(data)
    T = len(vocabulary)
    tables = {split: dataset.proteins(data, split) for split in ('development', 'test')}
    ids = [r['protein_id'] for r in tables['test']]
    inputs = scorer.prepared_inputs(scorer.data_checksums(data), vocabulary, tables, 'test', scorer_labels(labels), None,
                                    10_000, CACHE)
    if len(inputs['test_rows']) != len(ids):
        raise ValueError('every test protein must be scored (the label version drops proteins)')
    propagate = Propagation(parents)
    n = len(ids)
    y = torch.zeros(n * T, dtype=torch.bool, device=dev)
    y[torch.as_tensor(inputs['test'], device=dev)] = True
    y = y.view(n, T)
    left_out = torch.zeros(n * T, dtype=torch.bool, device=dev)
    left_out[torch.as_tensor(inputs['test_masked'], device=dev)] = True
    left_out = left_out.view(n, T)
    keep = ~left_out  # the pairs that count
    clusters, draws = inputs['clusters'], inputs['draws']
    bench = ev.Benchmark(y, ia, clusters, draws, device=dev, mask=left_out if left_out.any() else None)  # micro-AUPRC: the scorer's own code
    R, K = draws.shape

    # ---- macro set-up: positives sorted by term; per-cluster indicator; draws
    flat = torch.nonzero(y.reshape(-1)).squeeze(1)
    rows, terms = flat // T, flat % T
    order = torch.argsort(terms * n + rows)
    rows, terms = rows[order], terms[order]
    term_ids, per_term = torch.unique_consecutive(terms, return_counts=True)
    ends = torch.cumsum(per_term, 0) - 1
    starts = ends - per_term + 1
    M, NT = len(rows), len(term_ids)
    cl = torch.as_tensor(clusters, device=dev).long()
    onehot = torch.zeros(K, n, device=dev)
    onehot[cl, torch.arange(n, device=dev)] = 1
    W = torch.as_tensor(draws, dtype=torch.float32, device=dev)
    pos_cluster = cl[rows]
    say(f'test ({labels} labels): {n} proteins, {M} positive pairs on {NT} terms (the terms averaged), {int(left_out.sum())} left-out '
        f'pairs, {K} clusters, {R} replicates')

    def segment_sums(x):  # (r x M), positives sorted by term -> (r x NT) sums per term
        c = torch.cumsum(x, 1)
        return c[:, ends] - torch.where(starts > 0, c[:, (starts - 1).clamp_min(0)], torch.zeros_like(c[:, ends]))

    def macro(s, chunk_pairs=8192, chunk_reps=500):
        """(point, per-term AP, cluster replicates, terms averaged per replicate) of propagated scores s (n x T, float32 on the
        device)."""
        v = s[rows, terms]
        N = torch.empty(K, M, device=dev)
        TP = torch.empty(K, M, device=dev)
        for a in range(0, M, chunk_pairs):
            b = min(a + chunk_pairs, M)
            cols = terms[a:b]
            at_or_above = (s[:, cols] >= v[None, a:b]) & keep[:, cols]  # the positive's tie group and everything above it
            N[:, a:b] = onehot @ at_or_above.float()
            TP[:, a:b] = onehot @ (at_or_above & y[:, cols]).float()
        precision = TP.sum(0).double() / N.sum(0).double()
        ap = segment_sums(precision[None])[0] / per_term.double()
        reps, averaged = [], []
        for a in range(0, R, chunk_reps):
            w = W[a:a + chunk_reps]
            nw, tpw = (w @ N).double(), (w @ TP).double()
            wp = w[:, pos_cluster].double()
            contribution = torch.where(wp > 0, wp * tpw / nw.clamp_min(1), torch.zeros_like(wp))
            num, den = segment_sums(contribution), segment_sums(wp)
            has = den > 0
            reps.append((torch.where(has, num / den.clamp_min(1e-300), torch.zeros_like(num)).sum(1) / has.sum(1)))
            averaged.append(has.sum(1))
        return float(ap.mean()), ap.cpu().numpy(), torch.cat(reps).cpu().numpy(), torch.cat(averaged).cpu().numpy()

    def reference(s_np, y_np, keep_np, weights=None):
        """Mean over terms with a positive of the scorer's micro_auprc on each term's column (proteins repeated by weights; the
        pairs that do not count left out)."""
        if weights is not None:
            repeated = np.repeat(np.arange(len(s_np)), weights)
            s_np, y_np, keep_np = s_np[repeated], y_np[repeated], keep_np[repeated]
        values = [ev.micro_auprc(s_np[keep_np[:, t]][:, [t]], y_np[keep_np[:, t]][:, [t]]) for t in np.flatnonzero(y_np.any(0))]
        return float(np.mean(values)), len(values)

    def load(path):
        raw = open(path, 'rb').read()
        values = np.load(io.BytesIO(raw)).astype(np.float32)
        ids_file = Path(path).with_suffix('.ids.txt')
        if ids_file.exists():
            position = {p: i for i, p in enumerate(ids_file.read_text().split())}
            values = values[[position[p] for p in ids]]
        if values.shape != (n, T):
            raise ValueError(f'{path}: shape {values.shape}')
        return torch.as_tensor(values, device=dev), hashlib.sha256(raw).hexdigest()

    official = {}  # databank -> (test auprc values by (method, set, seed), replicates or None)
    for databank in {m['databank'] for m in models.values()}:
        values = {(r['method'], r['training_set'], r['seed']): float(r['value'])
                  for r in csv.DictReader(open(databank / 'metrics.tsv'), delimiter='\t')
                  if r['split'] == 'test' and r['metric'] == 'auprc'}
        z = np.load(databank / 'replicates.npz') if (databank / 'replicates.npz').exists() else None
        official[databank] = (values, z)
        say(f'databank {databank}: {len(values)} test AUPRC rows; replicates {"yes" if z is not None else "none (values checked only)"}')

    # ---- score
    seed_rows, point, reps, per_term_ap, worst = [], {}, {}, {}, dict(value=0.0, replicates=0.0)
    start, done, total = time.time(), 0, sum(len(m['files']) for m in models.values())
    for (method, training_set), m in models.items():
        values, z = official[m['databank']]
        per_seed = []
        for seed, path in m['files'].items():
            s, sha = load(path)
            propagate(s, inplace=True)
            micro = bench.score(s)
            micro_value, micro_reps = micro['metrics']['auprc'], micro['replicates'][:, 0]
            d = abs(micro_value - values[method, training_set, seed])
            worst['value'] = max(worst['value'], d)
            if z is not None:
                worst['replicates'] = max(worst['replicates'],
                                          float(np.abs(micro_reps - z[f'{method}/{training_set}/{seed}'][:, 0]).max()))
            if max(worst.values()) > 1e-9:
                raise ValueError(f'{method}/{training_set}/{seed}: micro-AUPRC {micro_value} does not reproduce the databank')
            value, ap, replicates, averaged = macro(s)
            if not done:  # the macro checks, on the first run
                say(f'cluster bootstrap: terms with a drawn positive per replicate, of {NT}: mean {averaged.mean():.0f}, '
                    f'2.5-97.5% {np.percentile(averaged, 2.5):.0f}-{np.percentile(averaged, 97.5):.0f}')
                s_np, y_np, keep_np = s.cpu().numpy(), y.cpu().numpy(), keep.cpu().numpy()
                ref, count = reference(s_np, y_np, keep_np)
                say(f'check macro point {method}/{training_set}/{seed}: {value:.12f}, reference {ref:.12f} over {count} terms')
                if abs(ref - value) > 1e-9 or count != NT:
                    raise ValueError('the macro point value does not reproduce the reference')
                for r in (0, 1):
                    ref, count = reference(s_np, y_np, keep_np, draws[r][clusters])
                    say(f'check macro replicate {r}: {replicates[r]:.12f}, reference {ref:.12f} over {count} terms')
                    if abs(ref - replicates[r]) > 1e-9:
                        raise ValueError('a macro replicate does not reproduce the reference')
            per_seed.append((micro_value, micro_reps, value, replicates, ap))
            seed_rows.append(dict(method=method, training_set=training_set, seed=seed, micro_auprc=micro_value,
                                  macro_auprc=value, path=path, sha256=sha))
            done += 1
            rate = done / (time.time() - start)
            print(f'{done}/{total} runs, {rate:.2f}/s, about {(total - done) / rate:.0f} s left: {method} {training_set} {seed} '
                  f'micro {micro_value:.4f} macro {value:.4f}', flush=True)
        key = (method, training_set)
        point['micro', key] = float(np.mean([p[0] for p in per_seed]))
        point['macro', key] = point['macro_cluster', key] = float(np.mean([p[2] for p in per_seed]))
        reps['micro', key] = np.mean([p[1] for p in per_seed], 0)
        reps['macro_cluster', key] = np.mean([p[3] for p in per_seed], 0)
        per_term_ap[key] = np.mean([p[4] for p in per_seed], 0)
        mean_reps = z[f'{method}/{training_set}/mean'][:, 0] if z is not None else None
        if mean_reps is not None:
            worst['replicates'] = max(worst['replicates'], float(np.abs(reps['micro', key] - mean_reps).max()))
        d = abs(point['micro', key] - values[method, training_set, 'mean'])
        worst['value'] = max(worst['value'], d)
        if max(worst.values()) > 1e-9:
            raise ValueError(f'{key}: the seed mean does not reproduce the databank')
    say(f'check micro-AUPRC against the databanks: every run and seed mean, largest difference {worst["value"]:.1e} (values), '
        f'{worst["replicates"]:.1e} (replicates)')

    # ---- the bootstrap over terms (the macro intervals): the same term draws for every model
    rng = np.random.default_rng(ev.BOOTSTRAP_SEED)
    counts = np.stack([np.bincount(rng.integers(0, NT, NT), minlength=NT) for _ in range(R)]).astype(np.float64)
    for key, ap in per_term_ap.items():
        if abs(ap.mean() - point['macro', key]) > 1e-12:
            raise ValueError(f'{key}: the seed-mean per-term values do not average to the macro point value')
        reps['macro', key] = counts @ ap / NT

    # ---- tables
    with open(out / 'seeds.tsv', 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(seed_rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(seed_rows)
    interval = lambda kind, key: tuple(float(x) for x in np.percentile(reps[kind, key], [2.5, 97.5]))  # noqa: E731
    with open(out / 'models.tsv', 'w') as f:
        f.write('method\ttraining_set\truns\tmicro_auprc\tmicro_ci_low\tmicro_ci_high\tmacro_auprc\tmacro_ci_low\tmacro_ci_high\t'
                'macro_run_sd\tmacro_terms\tmacro_cluster_ci_low\tmacro_cluster_ci_high\tmacro_cluster_replicate_mean\n')
        for key in models:
            runs = [r['macro_auprc'] for r in seed_rows if (r['method'], r['training_set']) == key]
            sd = f'{float(np.std(runs, ddof=1))!r}' if len(runs) > 1 else ''
            (lo_i, hi_i), (lo_a, hi_a), (lo_c, hi_c) = (interval(k, key) for k in ('micro', 'macro', 'macro_cluster'))
            f.write(f"{key[0]}\t{key[1]}\t{len(runs)}\t{point['micro', key]!r}\t{lo_i!r}\t{hi_i!r}\t{point['macro', key]!r}\t"
                    f"{lo_a!r}\t{hi_a!r}\t{sd}\t{NT}\t{lo_c!r}\t{hi_c!r}\t{float(reps['macro_cluster', key].mean())!r}\n")
    np.savez(out / 'replicates.npz', **{f'{kind}/{k[0]}/{k[1]}': v for (kind, k), v in reps.items()})
    np.savez(out / 'per_term.npz', terms=np.array([vocabulary[t] for t in term_ids.tolist()]),
             positives=per_term.cpu().numpy(), **{f'{k[0]}/{k[1]}': v for k, v in per_term_ap.items()})

    PUBLISHED = ('deepgose', 'transfew', 'dpfunc')
    ORIGINAL = {'deepgose': 'original_deepgose', 'transfew': 'original_transfew', 'dpfunc': 'original_dpfunc'}
    COMPARISONS = [('Published, original to HQ', (m, 'hq'), (m, 'original'), '') for m in PUBLISHED] + \
                  [('Published, HQ to ArchipelaGO', (m, 'archipelago'), (m, 'hq'), '') for m in PUBLISHED] + \
                  [('Published, original to ArchipelaGO', (m, 'archipelago'), (m, 'original'), 'yes') for m in PUBLISHED] + \
                  [('Finch, ArchipelaGO minus HQ + random', ('finch', 'archipelago'), ('finch', 'hq_random'), 'yes'),
                   ('Finch, ArchipelaGO minus HQ', ('finch', 'archipelago'), ('finch', 'hq'), '')] + \
                  [('Finch, ArchipelaGO minus an original set', ('finch', 'archipelago'), ('finch', ORIGINAL[m]), '')
                   for m in PUBLISHED] + \
                  [('Finch minus a published method on ArchipelaGO', ('finch', 'archipelago'), (m, 'archipelago'), 'yes')
                   for m in PUBLISHED]

    def contrast(kind, a, b):
        lo, hi = np.percentile(reps[kind, a] - reps[kind, b], [2.5, 97.5])
        return point[kind, a] - point[kind, b], float(lo), float(hi)

    header = ('comparison\ttable2\tmetric\tbootstrap\tmethod\ttraining_set\tminus_method\tminus_training_set\tvalue\tminus_value\t'
              'difference\tci_low\tci_high\n')
    KINDS = {'micro': ('micro_auprc', 'clusters'), 'macro': ('macro_auprc', 'terms'),
             'macro_cluster': ('macro_auprc', 'clusters (diagnostic, not used)')}
    with open(out / 'contrasts.tsv', 'w') as f:
        f.write(header)
        for kind, (metric, bootstrap) in KINDS.items():
            rows_out = list(COMPARISONS)
            on_archipelago = {m: point[kind, (m, 'archipelago')] for m in PUBLISHED}
            best, worst_m = sorted(on_archipelago, key=on_archipelago.get, reverse=True)[::2]
            rows_out.append(('Spread of the published methods on ArchipelaGO (best minus worst)', (best, 'archipelago'),
                             (worst_m, 'archipelago'), ''))
            for name, a, b, t2 in rows_out:
                d, lo, hi = contrast(kind, a, b)
                f.write(f'{name}\t{t2}\t{metric}\t{bootstrap}\t{a[0]}\t{a[1]}\t{b[0]}\t{b[1]}\t{point[kind, a]!r}\t'
                        f'{point[kind, b]!r}\t{d!r}\t{lo!r}\t{hi!r}\n')
    with open(out / 'contrasts_all.tsv', 'w') as f:
        f.write(header)
        keys = list(models)
        for kind in ('micro', 'macro'):
            metric, bootstrap = KINDS[kind]
            for a, b in itertools.permutations(keys, 2):
                if a[0] == b[0] or a[1] == b[1]:
                    d, lo, hi = contrast(kind, a, b)
                    f.write(f'\t\t{metric}\t{bootstrap}\t{a[0]}\t{a[1]}\t{b[0]}\t{b[1]}\t{point[kind, a]!r}\t{point[kind, b]!r}\t'
                            f'{d!r}\t{lo!r}\t{hi!r}\n')
    say(f'wrote {out} in {time.time() - start:.0f} s')
    log.close()


if __name__ == '__main__':
    main()
