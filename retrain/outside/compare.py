"""Finch and the published methods on each method's own benchmark: one table.

    python -m retrain.outside.compare {deepgose,transfew,dpfunc} --data DIR --replication REP
        -> DIR/outside/compare/{check,table,contrasts,relatives,coverage}.tsv

DIR is the method's outside root (retrain.outside.build; Finch combined in runs/finch/), REP its replication root
(retrain.replicate; the method's own retrain in runs/<method>/original/, its authors' model in replicate/released/).
Models: Finch arm A (original_<method>, seeds 9301-9303) and arm B (archipelago_screened); the method's own-data retrain
(DeepGO-SE's six-world ensemble, TransFew's last epoch seeds 0-2, DPFunc seeds 0-2); the authors' model (DeepGO-SE's
released models, TransFew's released prediction files, DPFunc's released checkpoints on our inputs).

Metrics, from saved predictions only:
  Fmax         the paper's protocol (DeepGO-SE compute_metrics, TransFew's CAFA-evaluator 1.2.1 call, DPFunc's
               new_compute_performance_deepgoplus), recomputed from per-protein counts taken from each paper's own scoring
               structures; check.tsv compares them with every official result available (retrain.outside.score,
               retrain.replicate.score; their AUPR when written with --paper-aupr).
  micro-AUPRC  ours: archipelago.evaluate (the paper's evaluator), average precision over all protein-term pairs pooled, on
               the propagated prediction and truth matrices each paper's code scores, restricted to its vocabulary without
               the MF root: DeepGO-SE predictions propagated over go.obo with every relation (the ensemble and the released
               models as saved, already propagated) and truth prop_annotations; TransFew predict.py's terms and > 0.01
               cutoff at full precision (the authors' files at their three decimals), CAFA-evaluator's max-propagation and
               truth; DPFunc its max-propagated predictions and closed truth. Proteins without a positive among the
               evaluated terms are left out of micro-AUPRC (archipelago.evaluate.Benchmark needs one).
Subsets: all; no_identical_training (no identical sequence in the method's training split); no_relative (no relative in it
under the paper's screen). Terms: all, or shared (training positives in both Finch arms).
Intervals: paired bootstrap, 10,000 replicates (archipelago.evaluate.cluster_draws, its fixed seed) over unique test
sequences (records with one sequence move together; none of the three benchmarks ships homology clusters of its test
set), the same draws for every model and both metrics, each model's statistic the mean over its seeds in every replicate;
2.5th and 97.5th percentiles of B - A, A - method and B - method.
relatives.tsv: test proteins with an identical sequence, a DIAMOND-only relative or none in arm A's training split.
coverage.tsv: at the highest threshold of the paper's grid, protein-centric recall and precision, the share of predicted
pairs (score > 0.01) at or above it, and the paper's trapezoid AUPR with and without the area from recall 0.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from retrain.common import add_argument, data_root, log
from retrain.outside.build import ARM_A, ARM_A_SCREENED, ARM_B
from retrain.replicate import score as replicate

ROOT_TERM = 'GO:0003674'
FINCH = (9301, 9302, 9303)
OWN = {'deepgose': ((0,), 'scores.npy'), 'transfew': ((0, 1, 2), 'scores.npy'), 'dpfunc': ((0, 1, 2), 'scores_raw.npy')}
RELEASED = {'deepgose': 'replicate/released/test_scores.npy', 'dpfunc': 'replicate/released/test_scores_raw.npy'}
CONTRASTS = (('finch_B', 'finch_A'), ('finch_A', 'method'), ('finch_B', 'method'),
             ('finch_B', 'finch_A_screened'), ('finch_A', 'finch_A_screened'))  # the last two where A-screened was trained


def read_tsv(path):
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


class Bench:
    """A method's test records, subsets and model scores (rows = the root's test dataset order, columns = its vocabulary)."""

    def __init__(self, method, root, replication):
        from archipelago.io import read_table
        self.method = method
        self.path, self.replication = Path(root), Path(replication)
        self.vocabulary = json.loads((self.path / 'datasets/ontology/vocabulary.json').read_text())
        self.shared = set(json.loads((self.path / 'outside/shared_terms.json').read_text()))
        self.ids = [r['seq_id'] for r in read_table(self.path / 'datasets/test/proteins.tsv')]
        rows = read_tsv(self.replication / 'replicate' / 'records.tsv')
        train = {r['seq_id'] for r in rows if r['split'] == 'train'}
        self.records = [r for r in rows if r['split'] == 'test']
        relative = {r['seq_id']: r['relative_in_original'] for r in read_tsv(self.path / 'outside/relatives.tsv')}
        for r in self.records:
            r['relative'] = relative[r['seq_id']]
            r['identical_training'] = r['seq_id'] in train
        self.models = {'finch_A': [('finch', ARM_A.format(method), s) for s in FINCH],
                       'finch_B': [('finch', ARM_B, s) for s in FINCH],
                       'method': [('own', None, s) for s in OWN[method][0]],
                       'released': [('released', None, 0)]}
        if (self.path / 'runs/finch' / ARM_A_SCREENED.format(method)).exists():  # arm A screened (TransFew, DPFunc)
            self.models['finch_A_screened'] = [('finch', ARM_A_SCREENED.format(method), s) for s in FINCH]

    def raw(self, kind, name, seed):
        """{seq_id: test score row} of one model run."""
        if kind == 'finch':
            path = self.path / 'runs/finch' / name / f'seed{seed}/test_scores.npy'
        elif kind == 'own':
            path = self.replication / 'runs' / self.method / 'original' / f'seed{seed}' / f'test_{OWN[self.method][1]}'
        else:
            path = self.replication / RELEASED[self.method]
        x = np.load(path)
        if x.shape != (len(self.ids), len(self.vocabulary)):
            raise ValueError(f'{path}: {x.shape}')
        return dict(zip(self.ids, x))


class Counts:
    """Per-protein counts of one model run under the paper's Fmax protocol: tp, npred (rows x thresholds), ntrue, and the
    rows the protocol counts (include)."""

    def __init__(self, protocol, tp, npred, ntrue, include):
        self.protocol, self.tp, self.npred = protocol, tp.astype(np.float64), npred.astype(np.float64)
        self.ntrue, self.include = ntrue.astype(np.float64), include

    def curve(self, w):
        """(precision, recall, valid points) with weights w (rows,) or (replicates, rows), the paper's averaging."""
        w = np.atleast_2d(w) * self.include
        has = self.npred > 0
        precision = np.divide(self.tp, self.npred, out=np.zeros_like(self.tp), where=has)
        recall = np.divide(self.tp, self.ntrue[:, None], out=np.zeros_like(self.tp), where=self.ntrue[:, None] > 0)
        covered = w @ has
        p = np.divide(w @ precision, covered, out=np.zeros_like(covered), where=covered > 0)
        r = (w @ recall) / w.sum(1, keepdims=True)
        valid = covered > 0 if self.protocol == 'transfew' else np.ones_like(covered, dtype=bool)  # cov > 0 rows only
        return p, r, valid

    def fmax(self, w):
        p, r, valid = self.curve(w)
        f = np.divide(2 * p * r, p + r, out=np.zeros_like(p), where=(p + r) > 0) * valid
        return f.max(1)

    def aupr(self, w):
        """The paper's protein-centric trapezoid AUPR (one weight vector) and the same with the area from recall 0."""
        p, r, valid = (a[0] for a in self.curve(w))
        p, r = p[valid], r[valid]
        order = np.argsort(r)  # the papers' own (default) sort of the points in threshold order
        area = float(np.trapz(p[order], r[order]))
        return area, area + float(p[order][0] * r[order][0])


# ---------------------------------------------------------------------------------------------------------- per paper
def deepgose(b, with_counts=True):
    """{term set: {(model, seed): Counts}}, micro (truth, {(model, seed): scores}, columns), rows, thresholds.
    with_counts=False: the micro-AUPRC matrices only (retrain.outside.masked)."""
    from retrain.deepgose.train import OfficialPropagation
    _, _, go, terms, _, frames = replicate.deepgose_protocol(b.replication)
    test_df = frames['test']
    log('deepgose test_data.pkl columns', list(test_df.columns))
    terms = list(terms)
    propagate = OfficialPropagation(b.replication / 'external' / 'deepgose' / 'data' / 'go.obo', terms)
    go_set = go.get_namespace_terms('molecular_function') - {ROOT_TERM}
    by_record = {r['record']: r for r in b.records}
    df = test_df[test_df['proteins'].isin(by_record)].reset_index(drop=True)
    rows = [by_record[p] for p in df['proteins']]
    index = {t: j for j, t in enumerate(terms)}
    in_go = np.array([t in go_set for t in terms])
    scores = {}
    for model, runs in b.models.items():
        for kind, name, seed in runs:
            raw = b.raw(kind, name, seed)
            x = np.stack([raw[r['seq_id']] for r in rows]).astype(np.float32)
            scores[model, seed] = propagate(x) if kind == 'finch' else x  # the ensemble and released models: as saved
    counts = {}
    for term_set in ('all', 'shared') if with_counts else ():
        keep = None if term_set == 'all' else b.shared
        labels = [{t for t in a if t in go_set and (keep is None or t in keep)} for a in df['prop_annotations']]
        cols = np.array([j for j, t in enumerate(terms) if keep is None or t in keep])
        lab = np.zeros((len(rows), len(cols)), bool)
        for i, a in enumerate(labels):
            for t in a:
                if t in index and (keep is None or t in keep):
                    lab[i, np.searchsorted(cols, index[t])] = True
        ntrue = np.array([len(a) for a in labels])
        g = in_go[cols]
        counts[term_set] = {}
        for key, x in scores.items():
            xs = x[:, cols]
            tp, npred = np.zeros((len(rows), 101), np.int64), np.zeros((len(rows), 101), np.int64)
            for t in range(101):
                p = (xs >= np.float32(t / 100.0)) & g
                npred[:, t], tp[:, t] = p.sum(1), (p & lab).sum(1)
            counts[term_set][key] = Counts('deepgose', tp, npred, ntrue, ntrue > 0)
    micro_cols = [j for j, t in enumerate(terms) if t in go_set]
    truth = np.zeros((len(rows), len(terms)), bool)
    for i, a in enumerate(df['prop_annotations']):
        for t in a:
            if t in index:
                truth[i, index[t]] = True
    micro = (truth[:, micro_cols], {k: v[:, micro_cols] for k, v in scores.items()}, [terms[j] for j in micro_cols])
    return {'test': (rows, counts, micro, 100)}


def transfew(b, with_counts=True):
    import cafaeval.parser as cp
    from retrain.transfew.original import index_pickle
    TFEW = replicate.tfew(b.replication)
    TFEW_GT = replicate.transfew_groundtruth(b.replication)
    TFEW_AUTHORS = {test_set: TFEW / 'evaluation' / 'predictions' / folder / 'swissprot_mf' / 'transfew.tsv'
                    for test_set, folder in (('Test_all', 'full'), ('Test_novel', 'seq_ID_30'))}
    transfew_file = replicate.transfew_file
    sorted_terms = index_pickle(TFEW, 'sorted_terms')
    full = list(index_pickle(TFEW, 'term_indicies')[0])
    ontology = cp.obo_parser(str(TFEW / 'obo/go-basic.obo'), ('is_a', 'part_of'), str(TFEW / 'test/output_t1_t2/ia_prop.txt'), True)
    mf = ontology['molecular_function']
    tau = np.arange(0.01, 1, 0.01)
    work = b.path / 'outside' / 'compare' / 'work'
    work.mkdir(parents=True, exist_ok=True)
    by_accession = {r['accession']: r for r in b.records}
    out = {}
    for test_set, gt_all in TFEW_GT.items():
        listed = {line.split('\t')[0] for line in open(gt_all)}
        proteins = [(r['accession'], r['seq_id']) for r in b.records if r['accession'] in listed]
        counts, rows = {}, None
        for term_set in ('all', 'shared') if with_counts else ('all',):
            keep = None if term_set == 'all' else b.shared
            cols = [j for j in full if keep is None or sorted_terms[j] in keep]
            terms = [sorted_terms[j] for j in cols]
            lines = [line for line in open(gt_all) if keep is None or line.rstrip('\n').split('\t')[1] in keep]
            if keep is not None:  # as evaluate.py: proteins left with only the root are left out
                counted = {line.split('\t')[0] for line in lines if line.rstrip('\n').split('\t')[1] != ROOT_TERM}
                lines = [line for line in lines if line.split('\t')[0] in counted]
            gt_file = work / f'gt_{test_set}_{term_set}.tsv'
            gt_file.write_text(''.join(lines))
            gt = cp.gt_parser(str(gt_file), ontology)['molecular_function']
            order = sorted(gt.ids, key=gt.ids.get)
            g = gt.matrix[:, mf.toi]
            counts[term_set] = {}
            for model, runs in b.models.items() if with_counts else ():
                for kind, name, seed in runs:
                    path = work / f'{model}_{seed}_{test_set}_{term_set}.tsv'
                    if kind == 'released':
                        src = TFEW_AUTHORS[test_set]
                        path.write_text(''.join(line for line in open(src) if keep is None or
                                                (line.count('\t') >= 2 and line.split('\t')[1] in keep)))
                    else:
                        transfew_file(path, b.raw(kind, name, seed), proteins, terms, cols)
                    p = cp.pred_parser(str(path), ontology, {'molecular_function': gt}, 'max')['molecular_function'].matrix[:, mf.toi]
                    tp, npred = np.zeros((len(order), len(tau)), np.int64), np.zeros((len(order), len(tau)), np.int64)
                    for k, t in enumerate(tau):
                        hit = p >= t
                        npred[:, k], tp[:, k] = hit.sum(1), (hit & g).sum(1)
                    counts[term_set][model, seed] = Counts('transfew', tp, npred, g.sum(1), np.ones(len(order), bool))
            if term_set == 'all':
                rows = [by_accession[a] for a in order]
                gt_rows = order
            else:  # the shared rows are a subset of the full rows: align them to the full row order
                pos = {a: i for i, a in enumerate(order)}
                for key, c in counts['shared'].items():
                    idx = [pos.get(a) for a in gt_rows]
                    take = lambda m, fill: np.stack([m[i] if i is not None else fill for i in idx])  # noqa: E731
                    counts['shared'][key] = Counts('transfew', take(c.tp, np.zeros(len(tau))), take(c.npred, np.zeros(len(tau))),
                                                   np.array([c.ntrue[i] if i is not None else 0 for i in idx]),
                                                   np.array([i is not None for i in idx]))
        # micro-AUPRC: the full truth; predictions at full precision (predict.py's terms, > 0.01), max-propagated
        gt = cp.gt_parser(str(work / f'gt_{test_set}_all.tsv'), ontology)['molecular_function']
        position = {t: mf.terms_dict[t]['index'] for t in sorted_terms if t in mf.terms_dict}
        micro_terms = [t for t in b.vocabulary if t != ROOT_TERM and t in position]
        micro_cols = [position[t] for t in micro_terms]
        written = np.array([j for j in full if sorted_terms[j] in position])  # predict.py's terms (all in this ontology)
        if len(written) != len(full):
            raise ValueError('a predict.py term is missing from TransFew\'s go-basic.obo')
        target = np.array([position[sorted_terms[j]] for j in written])
        scores = {}
        for model, runs in b.models.items():
            for kind, name, seed in runs:
                if kind == 'released':
                    m = cp.pred_parser(str(TFEW_AUTHORS[test_set]), ontology, {'molecular_function': gt}, 'max')['molecular_function'].matrix
                else:
                    raw = b.raw(kind, name, seed)
                    m = np.zeros(gt.matrix.shape, dtype='float')
                    for acc, i in gt.ids.items():
                        values = raw[by_accession[acc]['seq_id']][written].astype(float)
                        keep_ = values > 0.01
                        m[i, target[keep_]] = values[keep_]
                    cp.propagate(m, mf, mf.order, mode='max')
                scores[model, seed] = m[:, micro_cols].astype(np.float32)
        out[test_set] = (rows, counts, (gt.matrix[:, micro_cols], scores, micro_terms), len(tau) - 1)
    return out


def dpfunc(b, with_counts=True):
    import pandas as pd
    evaluation = replicate.dpfunc_evaluation(b.replication)
    DPFUNC_GO = b.replication / 'datasets' / 'ontology' / 'go-basic.obo'
    captured, fmax_function, namespace_terms = {}, evaluation.fmax, evaluation.Ontology.get_namespace_terms

    def capture_fmax(targets, scores):
        captured['targets'], captured['scores'] = np.asarray(targets), np.asarray(scores)
        return fmax_function(targets, scores)

    def capture_terms(self, namespace):
        captured['go_set'] = namespace_terms(self, namespace)
        return captured['go_set']
    evaluation.fmax = capture_fmax
    evaluation.Ontology.get_namespace_terms = capture_terms
    truth = replicate.dpfunc_truth(b.replication, 'test')  # mf_test2_go.txt within the output terms
    go = evaluation.Ontology(str(DPFUNC_GO), with_rels=True)
    column = {t: j for j, t in enumerate(b.vocabulary)}
    rows = list(b.records)
    raws = {(model, seed): b.raw(kind, name, seed) for model, runs in b.models.items() for kind, name, seed in runs}
    counts, micro = {}, None
    for term_set in ('all', 'shared') if with_counts else ('all',):
        keep = None if term_set == 'all' else b.shared
        include = np.array([keep is None or any(t != ROOT_TERM and go.has_term(t) for t in truth[r['record']] & keep)
                            for r in rows])
        chosen = [r for r, k in zip(rows, include) if k]
        names = [t for t in b.vocabulary if keep is None or t in keep]
        cols = [column[t] for t in names]
        counts[term_set] = {}
        for key, raw in raws.items():
            frame = pd.DataFrame(dict(protein_id=[r['record'] for r in chosen],
                                      gos=[truth[r['record']] if keep is None else truth[r['record']] & keep for r in chosen],
                                      predictions=[dict(zip(names, map(float, raw[r['seq_id']][cols]))) for r in chosen]))
            evaluation.new_compute_performance_deepgoplus(frame, str(DPFUNC_GO), 'mf')
            t, s = captured['targets'] > 0, captured['scores']
            tp, npred = np.zeros((len(rows), 101), np.int64), np.zeros((len(rows), 101), np.int64)
            sel = np.flatnonzero(include)
            for c in range(101):
                hit = s >= c / 100
                npred[sel, c], tp[sel, c] = hit.sum(1), (hit & t).sum(1)
            ntrue = np.zeros(len(rows))
            ntrue[sel] = t.sum(1)
            counts[term_set][key] = Counts('dpfunc', tp, npred, ntrue, include)
            if term_set == 'all':  # this call's column order (list(go_set), as the function made it)
                where = {x: j for j, x in enumerate(captured['go_set'])}
                if micro is None:
                    micro_terms = [x for x in b.vocabulary if x != ROOT_TERM and x in where]
                    micro = (t[:, [where[x] for x in micro_terms]], {}, micro_terms)
                cols_now = [where[x] for x in micro[2]]
                if not np.array_equal(t[:, cols_now], micro[0]):
                    raise ValueError('DPFunc truth differs between calls')
                micro[1][key] = np.clip(s[:, cols_now], 0, 1).astype(np.float32)
    return {'test': (rows, counts, micro, 100)}


# ------------------------------------------------------------------------------------------------------------- analysis
def official(b):
    """{(model, seed, test_set, subset, terms): (Fmax, AUPR or None)} of the official results available: retrain.replicate.score's
    (at 4 decimals), then retrain.outside.score's unrounded ones in their place where both exist."""
    import re
    out = {}
    paths = ((b.replication / 'replicate' / 'results' / 'scores.tsv', None),
             (b.path / 'outside' / 'results' / 'scores.tsv', {'original': 'finch_A', 'archipelago': 'finch_B'}))
    for path, model_of in paths:
        if not path.exists():
            continue
        for r in read_tsv(path):
            aupr = float(r['AUPR']) if r.get('AUPR') else None
            if model_of is None:
                seed = re.match(r'seed (\d+)', r['model'])
                key = (('released', 0) if r['model'].startswith(('released', 'authors')) else ('method', 0)
                       if r['model'].startswith('ensemble') else ('method', int(seed.group(1)))
                       if seed and 'lowest' not in r['model'] else None)
                if key:
                    out[(*key, r.get('test_set', 'test'), r['subset'], 'all')] = (float(r['Fmax']), aupr)
            else:
                out[model_of[r['arm']], int(r['seed']), r['test_set'], r['subset'], r['terms']] = (float(r['Fmax']), aupr)
    return out


def subsets(rows):
    return {'all': [True] * len(rows), 'no_identical_training': [not r['identical_training'] for r in rows],
            'no_relative': [r['relative'] == 'none' for r in rows]}


def main():
    from archipelago.evaluate import Benchmark, cluster_draws, micro_auprc
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('method', choices=['deepgose', 'transfew', 'dpfunc'])
    parser.add_argument('--replication', required=True, help="the method's replication root")
    add_argument(parser)
    a = parser.parse_args()
    b = Bench(a.method, data_root(a.data), a.replication)
    views = globals()[a.method](b)
    ref = official(b)
    out = b.path / 'outside' / 'compare'
    out.mkdir(parents=True, exist_ok=True)
    check, table, contrasts, relatives, coverage = [], [], [], [], []
    for test_set, (rows, counts, (truth, scores, micro_terms), top) in views.items():
        shared_cols = np.array([t in b.shared for t in micro_terms])
        seq = np.array([r['seq_id'] for r in rows])
        for term_set in ('all', 'shared'):
            T = truth if term_set == 'all' else truth[:, shared_cols]
            S = scores if term_set == 'all' else {k: v[:, shared_cols] for k, v in scores.items()}
            for subset, mask in subsets(rows).items():
                mask = np.array(mask)
                w = mask.astype(float)
                # the official numbers this reproduces
                for (model, seed), c in counts[term_set].items():
                    key = (model, seed, test_set, subset, term_set)
                    if key in ref:
                        f, (pa, _) = float(c.fmax(w)[0]), c.aupr(w)
                        check.append(dict(model=model, seed=seed, test_set=test_set, subset=subset, terms=term_set,
                                          Fmax=f, Fmax_official=ref[key][0], AUPR=pa, AUPR_official=ref[key][1],
                                          largest=max(abs(f - ref[key][0]), abs(pa - ref[key][1]) if ref[key][1] is not None else 0)))
                # bootstrap units: unique sequences of the subset
                rows_in = np.flatnonzero(mask)
                index, draws = cluster_draws(seq[rows_in])
                W = np.zeros((len(draws), len(rows)))
                W[:, rows_in] = draws[:, index]
                labelled = rows_in[T[rows_in].any(1)]
                present, relabel = np.unique(index[np.isin(rows_in, labelled)], return_inverse=True)
                bench = Benchmark(T[labelled], clusters=relabel, draws=draws[:, present])
                stats = {}
                for model in b.models:
                    keys = [k for k in counts[term_set] if k[0] == model]
                    fm = np.array([counts[term_set][k].fmax(w)[0] for k in keys])
                    fm_rep = np.mean([counts[term_set][k].fmax(W) for k in keys], axis=0)
                    res = [bench.score(S[k][labelled]) for k in keys]
                    au = np.array([r['metrics']['auprc'] for r in res])
                    au_rep = np.mean([r['replicates'][:, 0] for r in res], axis=0)
                    stats[model] = dict(Fmax=(fm, fm_rep), micro_AUPRC=(au, au_rep))
                    table.append(dict(test_set=test_set, subset=subset, terms=term_set, n_terms=int(T.shape[1]), model=model,
                                      runs=len(keys), proteins=int(mask.sum()), proteins_micro=len(labelled),
                                      sequences=int(draws.shape[1]), Fmax_mean=fm.mean(), Fmax_sd=fm.std(ddof=1) if len(fm) > 1 else 0.0,
                                      micro_AUPRC_mean=au.mean(), micro_AUPRC_sd=au.std(ddof=1) if len(au) > 1 else 0.0))
                    log(a.method, table[-1])
                for x, y in CONTRASTS:
                    if x not in stats or y not in stats:
                        continue
                    for metric in ('Fmax', 'micro_AUPRC'):
                        (px, rx), (py, ry) = stats[x][metric], stats[y][metric]
                        low, high = np.percentile(rx - ry, [2.5, 97.5])
                        contrasts.append(dict(test_set=test_set, subset=subset, terms=term_set, contrast=f'{x} - {y}',
                                              metric=metric, difference=px.mean() - py.mean(), ci_low=low, ci_high=high,
                                              units='unique test sequences', n_units=int(draws.shape[1]),
                                              replicates=len(draws)))
        # arm A's relatives (all terms, all proteins of the test set)
        for group in ('identical', 'diamond', 'none'):
            mask = np.array([r['relative'] == group for r in rows])
            if not mask.any():
                continue
            labelled = mask & truth.any(1)
            for model in b.models:
                keys = [k for k in counts['all'] if k[0] == model]
                fm = np.array([counts['all'][k].fmax(mask.astype(float))[0] for k in keys])
                au = np.array([micro_auprc(scores[k][labelled], truth[labelled]) for k in keys])
                relatives.append(dict(test_set=test_set, group=group, model=model, proteins=int(mask.sum()),
                                      Fmax_mean=fm.mean(), Fmax_sd=fm.std(ddof=1) if len(fm) > 1 else 0.0,
                                      micro_AUPRC_mean=au.mean(), micro_AUPRC_sd=au.std(ddof=1) if len(au) > 1 else 0.0))
        # curve coverage at the top of the paper's grid (all proteins, all terms)
        ones = np.ones(len(rows))
        for (model, seed), c in counts['all'].items():
            p, r, valid = (x[0] for x in c.curve(ones))
            area, with_rectangle = c.aupr(ones)
            x = scores[model, seed]
            predicted = x[x > 0.01]
            coverage.append(dict(test_set=test_set, model=model, seed=seed, top_threshold=(top + 1) / 100 if a.method == 'transfew' else top / 100,
                                 recall_at_top=r[top], precision_at_top=p[top],
                                 share_of_predictions_at_or_above_top=float(np.mean(predicted >= ((top + 1) / 100 if a.method == 'transfew' else top / 100))),
                                 paper_AUPR=area, paper_AUPR_with_area_from_recall_0=with_rectangle))
    for name, data in (('check', check), ('table', table), ('contrasts', contrasts), ('relatives', relatives), ('coverage', coverage)):
        with open(out / f'{name}.tsv', 'w', newline='') as f:
            if data:
                w = csv.DictWriter(f, fieldnames=list(data[0]), delimiter='\t')
                w.writeheader()
                w.writerows(data)
    worst = max((c['largest'] for c in check), default=float('nan'))
    log(a.method, f'check: {len(check)} official rows, largest difference {worst:.2e}; wrote {out}')
    if not worst <= 1e-4:
        raise ValueError(f'the recomputed Fmax / AUPR differ from the official results by {worst}')


if __name__ == '__main__':
    main()
