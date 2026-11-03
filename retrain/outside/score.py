"""Score Finch's two arms on a published method's own test set with that paper's own evaluation code.

    python -m retrain.outside.score {deepgose,transfew,dpfunc} --data DIR --replication REP [--seeds 9301 9302 9303]
        -> DIR/outside/results/{scores.tsv, summary.tsv}

The protocols are retrain.replicate.score's (REP holds the authors' files under external/), on Finch's raw scores
(retrain.finch.combine) of both arms:
  DeepGO-SE  each run's scores propagated as the official pipeline propagates a model's predictions (go.obo, is_a and every
             relationship), then deepgo2's compute_metrics: Fmax, Smin, AUPR, AUC
  TransFew   the released predict.py's output, CAFA-evaluator 1.2.1 on Test_all and Test_novel: Fmax, WFmax, AUPR, Smin
  DPFunc     DPFunc's new_compute_performance_deepgoplus: Fmax, AUPR
F1: the test F at the threshold of each run's Fmax on the method's validation set (the same code on the authors' validation
proteins and labels); Fmax itself is test-selected. Rows per arm (original: the method's own training split; archipelago:
ArchipelaGO screened), seed, test set, subset (all; no_relative: the test proteins without a relative in either training set,
outside/relatives.tsv) and terms (all: the full vocabulary; shared: the terms with training positives in both arms, their
number in n_terms; truth and predictions restricted to them). summary.tsv: the mean and SD over seeds per arm, and of the
per-seed difference archipelago - original.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from archipelago.io import read_table
from retrain.common import add_argument, data_root, load, log
from retrain.replicate import score as replicate

ARMS = {'original': 'original_{}', 'archipelago': 'archipelago_screened'}


class Benchmark:
    """The outside root's validation and test proteins, the authors' records of them, the subsets and the shared terms."""

    def __init__(self, method, root, replication, seeds):
        self.method, self.root, self.replication, self.seeds = method, Path(root), Path(replication), seeds
        self.ids = {split: [r['seq_id'] for r in read_table(self.root / 'datasets' / name / 'proteins.tsv')]
                    for split, name in (('valid', 'development'), ('test', 'test'))}
        rows = replicate.records(replication)
        self.records = {split: [r for r in rows if r['split'] == split and r['seq_id'] in set(self.ids[split])]
                        for split in ('valid', 'test')}
        relative = {r['seq_id']: r['relative_in_original'] for r in read_table(self.root / 'outside' / 'relatives.tsv')}
        self.subsets = {'all': {r['record'] for r in self.records['test']},
                        'no_relative': {r['record'] for r in self.records['test'] if relative[r['seq_id']] == 'none'}}
        self.vocabulary = json.loads((self.root / 'datasets/ontology/vocabulary.json').read_text())
        self.shared = set(json.loads((self.root / 'outside' / 'shared_terms.json').read_text()))

    def scores(self, arm, seed, split):
        x = np.load(self.root / 'runs' / 'finch' / ARMS[arm].format(self.method) / f'seed{seed}' /
                    f"{'dev' if split == 'valid' else 'test'}_scores.npy")
        if x.shape != (len(self.ids[split]), len(self.vocabulary)):
            raise ValueError(f'{arm} seed {seed} {split}: scores of shape {x.shape}')
        return dict(zip(self.ids[split], x))


def f_of(p, r):
    return 2 * p * r / (p + r) if p + r > 0 else 0.0


def deepgose(b, rows):
    from retrain.deepgose.train import OfficialPropagation
    compute_metrics, evaluate_annotations, go, terms, terms_dict, frames = replicate.deepgose_protocol(b.replication)
    propagate = OfficialPropagation(b.replication / 'external' / 'deepgose' / 'data' / 'go.obo', list(terms))
    go_set = go.get_namespace_terms('molecular_function') - {'GO:0003674'}

    def frame(split, keep, subset_terms):
        recs = {r['record']: r['seq_id'] for r in b.records[split] if r['record'] in keep}
        df = frames[split][frames[split]['proteins'].isin(recs)].reset_index(drop=True).copy()
        if subset_terms is not None:
            df['prop_annotations'] = [[t for t in a if t in subset_terms] for a in df['prop_annotations']]
        return df, [recs[p] for p in df['proteins']]

    for arm in ARMS:
        for seed in b.seeds:
            propagated = {s: dict(zip(b.ids[s], propagate(np.stack(list(b.scores(arm, seed, s).values())))))
                          for s in ('valid', 'test')}
            for term_set, subset_terms in (('all', None), ('shared', b.shared)):
                idx = [terms_dict[t] for t in terms if subset_terms is None or t in subset_terms]
                sub_terms = np.asarray(terms)[idx]
                sub_dict = {t: i for i, t in enumerate(sub_terms)}
                vdf, vids = frame('valid', {r['record'] for r in b.records['valid']}, subset_terms)
                threshold = compute_metrics(vdf, go, sub_dict, sub_terms, 'mf',
                                            np.stack([propagated['valid'][h][idx] for h in vids]))[2]
                for subset, keep in b.subsets.items():
                    df, ids = frame('test', keep, subset_terms)
                    pred = np.stack([propagated['test'][h][idx] for h in ids])
                    fmax, smin, tmax, _, _, auc, aupr, _, _ = compute_metrics(df, go, sub_dict, sub_terms, 'mf', pred)
                    labels = [set(x for x in a if x in go_set) for a in df['prop_annotations']]
                    chosen = [{sub_terms[j] for j in np.flatnonzero(p >= threshold)} & go_set for p in pred]
                    rows.append(dict(arm=arm, seed=seed, test_set='test', subset=subset, terms=term_set, n_terms=len(idx),
                                     proteins=len(df), Fmax=fmax, F1=evaluate_annotations(go, labels, chosen)[0],
                                     threshold_valid=threshold, threshold_test=tmax, Smin=smin, AUPR=aupr, AUC=auc))
                    log(rows[-1])


def transfew(b, rows):
    terms_all, columns_all = replicate.transfew_terms(b.replication)
    out = b.root / 'outside' / 'results' / 'work'
    out.mkdir(parents=True, exist_ok=True)
    valid = load(b.root, 'development')
    valid_terms = {h: [b.vocabulary[j] for j in valid.labels[i].indices] for i, h in enumerate(valid.seq_ids)}
    gt_valid = out / 'gt_valid.tsv'
    with open(gt_valid, 'w') as handle:  # the validation records' released (propagated) terms as a ground-truth file
        handle.writelines(f"{r['accession']}\t{t}\n" for r in b.records['valid'] for t in valid_terms[r['seq_id']])
    tests = {name: (gt, {line.split('\t')[0] for line in open(gt)})
             for name, gt in replicate.transfew_groundtruth(b.replication).items()}

    def restricted(gt, path, subset_terms, accessions=None):
        path.write_text(''.join(line for line in open(gt)
                                if (subset_terms is None or line.rstrip('\n').split('\t')[1] in subset_terms)
                                and (accessions is None or line.split('\t')[0] in accessions)))
        return path
    for term_set, subset_terms in (('all', None), ('shared', b.shared)):
        pairs = [(t, j) for t, j in zip(terms_all, columns_all) if subset_terms is None or t in subset_terms]
        terms, columns = [t for t, _ in pairs], [j for _, j in pairs]
        for arm in ARMS:
            for seed in b.seeds:
                path = out / f'{arm}_{seed}_valid_{term_set}.tsv'
                replicate.transfew_file(path, b.scores(arm, seed, 'valid'),
                                        [(r['accession'], r['seq_id']) for r in b.records['valid']], terms, columns)
                table = replicate.cafa_table(b.replication, path, restricted(gt_valid, out / f'gt_valid_{term_set}.tsv',
                                                                             subset_terms), out)
                threshold = float(table.loc[table['f'].idxmax(), 'tau'])
                scores = b.scores(arm, seed, 'test')
                for test_set, (gt, listed) in tests.items():
                    for subset, keep in b.subsets.items():
                        proteins = [(r['accession'], r['seq_id']) for r in b.records['test']
                                    if r['record'] in keep and r['accession'] in listed]
                        gt_path = restricted(gt, out / f'gt_{test_set}_{subset}_{term_set}.tsv', subset_terms,
                                             {p for p, _ in proteins})
                        path = out / f'{arm}_{seed}_{test_set}_{subset}_{term_set}.tsv'
                        replicate.transfew_file(path, scores, proteins, terms, columns)
                        table = replicate.cafa_table(b.replication, path, gt_path, out)
                        at = table.iloc[(table['tau'] - threshold).abs().argmin()]
                        rows.append(dict(arm=arm, seed=seed, test_set=test_set, subset=subset, terms=term_set,
                                         n_terms=len(columns), proteins=len(proteins), **replicate.cafa_metrics(table),
                                         F1=float(at['f']), threshold_valid=threshold,
                                         threshold_test=float(table.loc[table['f'].idxmax(), 'tau'])))
                        log(rows[-1])


def dpfunc(b, rows):
    import pandas as pd
    evaluation = replicate.dpfunc_evaluation(b.replication)
    captured, official_fmax = {}, evaluation.fmax

    def capturing(targets, scores):  # the official fmax, keeping its precision and recall at every threshold
        out = official_fmax(targets, scores)
        captured['p'], captured['r'] = out[2], out[3]
        return out
    evaluation.fmax = capturing
    truth = {**replicate.dpfunc_truth(b.replication, 'valid'), **replicate.dpfunc_truth(b.replication, 'test')}
    column = {t: j for j, t in enumerate(b.vocabulary)}
    obo = str(b.root / 'datasets/ontology/go-basic.obo')

    def score(split, keep, scores, subset_terms):
        # with a term subset, proteins without a true term in it are left out (the official fmax averages recall over
        # proteins and needs a true term for each; DeepGO-SE's code and CAFA-evaluator leave such proteins out too)
        recs = [r for r in b.records[split] if r['record'] in keep
                and (subset_terms is None or truth[r['record']] & subset_terms)]
        names = [t for t in b.vocabulary if subset_terms is None or t in subset_terms]
        cols = [column[t] for t in names]
        frame = pd.DataFrame(dict(
            protein_id=[r['record'] for r in recs],
            gos=[truth[r['record']] if subset_terms is None else truth[r['record']] & subset_terms for r in recs],
            predictions=[dict(zip(names, map(float, scores[r['seq_id']][cols]))) for r in recs]))
        fmax, aupr, t = evaluation.new_compute_performance_deepgoplus(frame, obo, 'mf')
        return fmax, aupr, t, list(captured['p']), list(captured['r']), len(recs), len(cols)

    for arm in ARMS:
        for seed in b.seeds:
            for term_set, subset_terms in (('all', None), ('shared', b.shared)):
                threshold = score('valid', {r['record'] for r in b.records['valid']}, b.scores(arm, seed, 'valid'),
                                  subset_terms)[2]
                for subset, keep in b.subsets.items():
                    fmax, aupr, t, p, r, n, k = score('test', keep, b.scores(arm, seed, 'test'), subset_terms)
                    cut = int(round(threshold * 100))
                    rows.append(dict(arm=arm, seed=seed, test_set='test', subset=subset, terms=term_set, n_terms=k,
                                     proteins=n, Fmax=float(fmax), F1=f_of(p[cut], r[cut]), threshold_valid=threshold,
                                     threshold_test=t, AUPR=float(aupr)))
                    log(rows[-1])


def summarize(rows):
    """Seed means and SDs per arm, and of the per-seed difference archipelago - original."""
    import pandas as pd
    df = pd.DataFrame(rows)
    keys = ['test_set', 'subset', 'terms', 'n_terms', 'proteins']
    metrics = [c for c in df.columns if c not in keys + ['arm', 'seed', 'threshold_valid', 'threshold_test']]
    out = []
    for key, group in df.groupby(keys, sort=False):
        for arm in ARMS:
            g = group[group['arm'] == arm]
            out.append(dict(zip(keys, key), arm=arm, **{f'{m}_mean': g[m].mean() for m in metrics},
                            **{f'{m}_sd': g[m].std() for m in metrics}))
        d = group[group['arm'] == 'archipelago'].set_index('seed')[metrics] - \
            group[group['arm'] == 'original'].set_index('seed')[metrics]
        out.append(dict(zip(keys, key), arm='archipelago - original', **{f'{m}_mean': d[m].mean() for m in metrics},
                        **{f'{m}_sd': d[m].std() for m in metrics}))
    return pd.DataFrame(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('method', choices=['deepgose', 'transfew', 'dpfunc'])
    parser.add_argument('--replication', required=True, help="the method's replication root")
    parser.add_argument('--seeds', nargs='+', type=int, default=[9301, 9302, 9303])
    add_argument(parser)
    a = parser.parse_args()
    import pandas as pd
    b = Benchmark(a.method, data_root(a.data), a.replication, a.seeds)
    rows = []
    globals()[a.method](b, rows)
    folder = b.root / 'outside' / 'results'
    folder.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(folder / 'scores.tsv', sep='\t', index=False)
    if len(a.seeds) > 1:
        summarize(rows).to_csv(folder / 'summary.tsv', sep='\t', index=False)
    log(f'wrote {folder}')


if __name__ == '__main__':
    main()
