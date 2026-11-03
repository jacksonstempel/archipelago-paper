"""Score a replication root's test predictions with the paper's own evaluation code, and write the test set as scored.

    python -m retrain.replicate.score {deepgose,transfew,dpfunc} [--testset] [--paper-aupr] --data DIR
        -> replicate/results/scores.tsv (and with --testset: replicate/testset/)

Reported: each paper's Fmax (and Smin, AUC, WFmax where it reports them) by its own protocol, and the paper's (our)
micro-AUPRC on the matrices that protocol scores (retrain.replicate.micro). With --paper-aupr also each paper's own
protein-centric trapezoid AUPR, the check that its protocol is reproduced (it gives the papers' tables).

DeepGO-SE  deepgo2's deepgo/metrics.compute_metrics (the function its evaluate_entailment.py calls; external/deepgose/deepgo2
           at commit 6c77c98) on the 2,845 test records, with the ontology and information content exactly as
           evaluate_entailment.py builds them (go.obo with relationships; IC from the train, validation and test
           prop_annotations), on the six-world ensemble (retrain.deepgose.train ensemble: each world propagated, then
           averaged, as evaluate_entailment.py -c avg -nm 6). Fmax, Smin, AUPR, AUC. (evaluate_entailment.py itself reads
           nextprot_predictions_* whatever -td says, so the metric function is called directly.)
TransFew   CAFA-evaluator 1.2.1 (pip install cafaeval==1.2.1) as the authors' evaluate.py calls it (TFewData/obo/go-basic.obo,
           -ia test/output_t1_t2/ia_prop.txt, -norm cafa, -th_step 0.01, propagation max) on Test_all and Test_novel; AUPR =
           trapezoidal area under the evaluator's (recall, precision) points; Smin = the IA-weighted S. Predictions written as
           the released predict.py writes them (the terms of term_indicies[0], scores > 0.01, three decimals), for each seed
           and both kept checkpoints (last epoch, the primary; lowest validation BCE). The authors' released prediction
           files are scored the same way; they give the paper's Table 1 exactly.
DPFunc     DPFunc/evaluation.new_compute_performance_deepgoplus (external/dpfunc/DPFunc at commit c5faf4df) with the root's
           GO file, on each seed's mean of the last three epochs' sigmoid scores before post-processing (test_scores_raw.npy;
           the function applies the post-processing itself), with the test labels as DPFunc_main.py makes them.
Each method is also scored on its test proteins without an identical sequence in the authors' training split. Released
models (retrain.replicate.released) are scored when their scores exist.
"""
import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

import numpy as np

from archipelago.io import read_fasta, read_table
from retrain.common import add_argument, data_root, log
from retrain.replicate import micro

PAPER_AUPR = False  # --paper-aupr: also each paper's own protein-centric trapezoid AUPR (the papers' protocol check)

PAPER = {
    'deepgose': dict(source='Kulmanov et al. 2024, Nat Mach Intell, Table 1 (MF, UniProtKB/Swiss-Prot test set)',
                     Fmax=0.554, Smin=11.681, AUPR=0.552, AUC=0.874),
    'transfew': dict(source='Boadu and Cheng 2024, Bioinformatics Advances, Table 1 (MF); its Table 2 gives AUPR 0.6871 for '
                            'the same model', Test_all=dict(Fmax=0.6655, WFmax=0.5860, AUPR=0.3633, Smin=6.3848),
                     Test_novel=dict(Fmax=0.6476, WFmax=0.5601, AUPR=0.3513, Smin=6.8034)),
    'dpfunc': dict(source='Wang et al. 2025, Nat Commun, Table 2 (MF, large-scale dataset; mean of five runs)',
                   Fmax=0.635, AUPR=0.658),
}


def records(root):
    with open(Path(root) / 'replicate' / 'records.tsv', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def test_records(root):
    """The authors' test records with identical_training = 1 if the same sequence is in their training split."""
    rows = records(root)
    train = {r['seq_id'] for r in rows if r['split'] == 'train'}
    return [dict(r, identical_training=int(r['seq_id'] in train)) for r in rows if r['split'] == 'test']


def scores_by_protein(root, path, dataset='test'):
    """{seq_id: score row} of a score file whose rows follow datasets/<dataset>/proteins.tsv."""
    ids = [r['seq_id'] for r in read_table(Path(root) / 'datasets' / dataset / 'proteins.tsv')]
    x = np.load(path)
    if x.shape[0] != len(ids):
        raise ValueError(f'{path}: {x.shape[0]} rows for {len(ids)} {dataset} proteins')
    return dict(zip(ids, x))


def write_tsv(path, header, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)


# ------------------------------------------------------------------------------------------------------------ DeepGO-SE
def deepgose_protocol(root):
    """(compute_metrics, evaluate_annotations, go, terms, terms_dict, {split: dataframe}) as evaluate_entailment.py."""
    import pandas as pd
    data = Path(root) / 'external' / 'deepgose' / 'data'
    sys.path.insert(0, str(Path(root) / 'external' / 'deepgose' / 'deepgo2'))
    from deepgo.metrics import compute_metrics, evaluate_annotations
    from deepgo.utils import Ontology
    terms = pd.read_pickle(data / 'mf' / 'terms.pkl')['gos'].values.flatten()
    go = Ontology(str(data / 'go.obo'), with_rels=True)
    frames = {s: pd.read_pickle(data / 'mf' / f'{s}_data.pkl') for s in ('train', 'valid', 'test')}
    go.calculate_ic([set(x) for s in ('train', 'valid', 'test') for x in frames[s]['prop_annotations']])
    return compute_metrics, evaluate_annotations, go, terms, {t: i for i, t in enumerate(terms)}, frames


def deepgose(root, out):
    compute_metrics, _, go, terms, terms_dict, frames = deepgose_protocol(root)
    tests = {r['record']: r for r in test_records(root)}
    test = frames['test'][frames['test']['proteins'].isin(tests)].reset_index(drop=True)
    folder = Path(root) / 'runs' / 'deepgose' / 'original' / 'seed0'
    sources = {'ensemble': folder / 'test_scores.npy',
               'released models 0, 1, 2, 5, 6, 8': Path(root) / 'replicate' / 'released' / 'test_scores.npy'}
    rows = []
    for label, path in sources.items():
        if not path.exists():
            log(f'{path} not there; skipped')
            continue
        scores = scores_by_protein(root, path)
        go_set = go.get_namespace_terms('molecular_function') - {micro.ROOT}
        for subset in ('all', 'no_identical_training'):
            keep = [subset == 'all' or not tests[p]['identical_training'] for p in test['proteins']]
            if subset != 'all' and all(keep):
                continue
            df = test[keep].reset_index(drop=True)
            preds = np.stack([scores[tests[p]['seq_id']] for p in df['proteins']]).astype(np.float32)
            fmax, smin, tmax, _, _, auc, aupr, _, _ = compute_metrics(df, go, terms_dict, terms, 'mf', preds)
            rows.append([label, subset, len(df), round(fmax, 4), round(smin, 3),
                         round(micro.deepgose(preds, df['prop_annotations'], list(terms), go_set), 4), round(auc, 4), tmax,
                         *([round(aupr, 4)] if PAPER_AUPR else [])])
            log(rows[-1])
    write_tsv(out / 'scores.tsv', ['model', 'subset', 'proteins', 'Fmax', 'Smin', 'micro_AUPRC', 'AUC', 'threshold',
                                   *(['AUPR'] if PAPER_AUPR else [])], rows)


# ------------------------------------------------------------------------------------------------------------- TransFew
def tfew(root):
    return Path(root) / 'external' / 'transfew' / 'TFewData'


def transfew_groundtruth(root):
    t = tfew(root) / 'test' / 'output_t1_t2' / 'groundtruths'
    return dict(Test_all=t / 'full' / 'mf_swissprot.tsv', Test_novel=t / 'seq_ID_30' / 'mf_swissprot.tsv')


def transfew_terms(root):
    """predict.py's output terms: sorted_terms at term_indicies[0] (every term with a training annotation but the root)."""
    from retrain.transfew.original import index_pickle
    sorted_terms = index_pickle(tfew(root), 'sorted_terms')
    full = list(index_pickle(tfew(root), 'term_indicies')[0])
    return [sorted_terms[j] for j in full], full


def transfew_file(path, scores, proteins, terms, columns):
    """The released predict.py's output: per protein, terms with score > 0.01, '%0.3f', highest first."""
    with open(path, 'w') as handle:
        for accession, h in proteins:
            row = scores[h][columns]
            order = np.argsort(-row, kind='stable')
            handle.writelines('%s\t%s\t%0.3f\n' % (accession, terms[j], row[j]) for j in order if row[j] > 0.01)


def cafa_table(root, pred_file, gt_file, work):
    """CAFA-evaluator 1.2.1 as TransFew's evaluate.py runs it: its molecular_function rows over every threshold."""
    from cafaeval.evaluation import cafa_eval
    folder = Path(work) / 'pred'
    shutil.rmtree(folder, ignore_errors=True)
    folder.mkdir(parents=True)
    shutil.copyfile(pred_file, folder / 'prediction.tsv')
    df, _ = cafa_eval(str(tfew(root) / 'obo' / 'go-basic.obo'), str(folder), str(gt_file),
                      ia=str(tfew(root) / 'test' / 'output_t1_t2' / 'ia_prop.txt'), norm='cafa', prop='max', th_step=0.01,
                      n_cpu=4)
    d = df.reset_index()
    return d[d['ns'] == 'molecular_function'].sort_values('tau').reset_index(drop=True)


def cafa_metrics(table):
    by_recall = table.sort_values('rc')
    return dict(Fmax=float(table['f'].max()), WFmax=float(table['f_w'].max()),
                AUPR=float(np.trapz(by_recall['pr'], by_recall['rc'])), Smin=float(table['s_w'].min()))


def transfew(root, out):
    terms, columns = transfew_terms(root)
    tests = test_records(root)
    matrices = micro.TransFew(tfew(root), json.loads((Path(root) / 'datasets/ontology/vocabulary.json').read_text()))
    rows = []
    for test_set, gt in transfew_groundtruth(root).items():
        listed = {line.split('\t')[0] for line in open(gt)}
        for subset in ('all', 'no_identical_training'):
            proteins = [(r['accession'], r['seq_id']) for r in tests
                        if r['accession'] in listed and (subset == 'all' or not r['identical_training'])]
            keep = {a for a, _ in proteins}
            gt_file = gt
            if keep != listed:
                gt_file = out / f'gt_{test_set}_{subset}.tsv'
                gt_file.write_text(''.join(line for line in open(gt) if line.split('\t')[0] in keep))
            sources = {}  # label: (prediction file as predict.py writes it, the scores at full precision or None)
            if subset == 'all':
                sources["authors' released predictions"] = (tfew(root) / 'evaluation' / 'predictions' / {
                    'Test_all': 'full', 'Test_novel': 'seq_ID_30'}[test_set] / 'swissprot_mf' / 'transfew.tsv', None)
            for seed in range(3):
                fit = Path(root) / 'runs' / 'transfew' / 'original' / f'seed{seed}'
                for checkpoint, name in (('last epoch (primary)', 'test_scores.npy'),
                                         ('lowest validation BCE', 'test_scores_bcebest.npy')):
                    if (fit / name).exists():
                        path = out / 'predictions' / f'seed{seed}_{name[:-4]}_{test_set}_{subset}.tsv'
                        path.parent.mkdir(parents=True, exist_ok=True)
                        scores = scores_by_protein(root, fit / name)
                        transfew_file(path, scores, proteins, terms, columns)
                        sources[f'seed {seed}, {checkpoint}'] = (path, scores)
            for label, (path, scores) in sources.items():
                m = cafa_metrics(cafa_table(root, path, gt_file, out / 'work'))
                m['micro_AUPRC'] = (matrices.score(gt_file, pred_file=path) if scores is None else
                                    matrices.score(gt_file, None, scores, proteins, terms, columns))
                rows.append([label, test_set, subset, len(proteins), *(round(m[k], 4) for k in (
                    'Fmax', 'WFmax', 'micro_AUPRC', 'Smin', *(['AUPR'] if PAPER_AUPR else [])))])
                log(rows[-1])
    write_tsv(out / 'scores.tsv', ['model', 'test_set', 'subset', 'proteins', 'Fmax', 'WFmax', 'micro_AUPRC', 'Smin',
                                   *(['AUPR'] if PAPER_AUPR else [])], rows)


# --------------------------------------------------------------------------------------------------------------- DPFunc
def dpfunc_evaluation(root):
    """DPFunc/evaluation.py; get_anchestors (a pure breadth-first search) memoised for speed only."""
    import functools
    sys.path.insert(0, str(Path(root) / 'external' / 'dpfunc' / 'DPFunc'))
    from DPFunc import evaluation
    evaluation.Ontology.get_anchestors = functools.lru_cache(maxsize=None)(evaluation.Ontology.get_anchestors)
    return evaluation


def dpfunc_truth(root, split='test'):
    """{entry name: terms}: mf_test2_go.txt (split test) or mf_test1_go.txt (valid) within the 6,086 output terms, as
    mlb.transform gives DPFunc_main.py's labels."""
    from retrain.dpfunc.original import records as dpfunc_records
    vocabulary = set(json.loads((Path(root) / 'datasets/ontology/vocabulary.json').read_text()))
    return {name: set(terms) & vocabulary for name, s, _, terms in dpfunc_records(root) if s == split}


def dpfunc_frame(tests, scores, vocabulary, truth):
    import pandas as pd
    return pd.DataFrame(dict(protein_id=[r['record'] for r in tests], gos=[truth[r['record']] for r in tests],
                             predictions=[dict(zip(vocabulary, map(float, scores[r['seq_id']]))) for r in tests]))


def dpfunc(root, out):
    evaluation = dpfunc_evaluation(root)
    vocabulary = json.loads((Path(root) / 'datasets/ontology/vocabulary.json').read_text())
    truth, tests = dpfunc_truth(root), test_records(root)
    sources = {f'seed {k}': Path(root) / 'runs' / 'dpfunc' / 'original' / f'seed{k}' / 'test_scores_raw.npy' for k in range(3)}
    sources['released checkpoints'] = Path(root) / 'replicate' / 'released' / 'test_scores_raw.npy'
    rows = []
    for label, path in sources.items():
        if not path.exists():
            log(f'{path} not there; skipped')
            continue
        scores = scores_by_protein(root, path)
        for subset in ('all', 'no_identical_training'):
            chosen = [r for r in tests if subset == 'all' or not r['identical_training']]
            (fmax, aupr, t), (micro_auprc,) = micro.dpfunc(evaluation, dpfunc_frame(chosen, scores, vocabulary, truth),
                                                           str(Path(root) / 'datasets/ontology/go-basic.obo'), vocabulary)
            rows.append([label, subset, len(chosen), round(float(fmax), 4), round(micro_auprc, 4), t,
                         *([round(float(aupr), 4)] if PAPER_AUPR else [])])
            log(rows[-1])
    write_tsv(out / 'scores.tsv', ['model', 'subset', 'proteins', 'Fmax', 'micro_AUPRC', 'threshold',
                                   *(['AUPR'] if PAPER_AUPR else [])], rows)


# ------------------------------------------------------------------------------------------------------------- test set
def testset(root, method):
    """replicate/testset/: proteins.tsv (record, accession, seq_id, identical_training[, novel]), sequences.fasta (header:
    the authors' record id and the seq_id), labels.tsv (record, go_id: the truth exactly as handed to the paper's scoring
    code), vocabulary.txt, and the scoring files (GO file; TransFew's IA file and Test_novel ground truth)."""
    folder = Path(root) / 'replicate' / 'testset'
    folder.mkdir(parents=True, exist_ok=True)
    tests = test_records(root)
    sequences = read_fasta(Path(root) / 'datasets' / 'test' / 'sequences.fasta')
    extra = ['novel'] if method == 'transfew' else []
    write_tsv(folder / 'proteins.tsv', ['record', 'accession', 'seq_id', 'identical_training', *extra],
              [(r['record'], r['accession'], r['seq_id'], r['identical_training'], *[r[c] for c in extra]) for r in tests])
    with open(folder / 'sequences.fasta', 'w') as handle:
        handle.writelines(f">{r['record']} {r['seq_id']}\n{sequences[r['seq_id']]}\n" for r in tests)
    vocabulary = json.loads((Path(root) / 'datasets/ontology/vocabulary.json').read_text())
    (folder / 'vocabulary.txt').write_text(''.join(t + '\n' for t in vocabulary))
    if method == 'deepgose':  # compute_metrics: prop_annotations within the MF namespace of go.obo, root removed
        _, _, go, _, _, frames = deepgose_protocol(root)
        go_set = go.get_namespace_terms('molecular_function') - {'GO:0003674'}
        keep = {r['record'] for r in tests}
        labels = [(p, t) for p, a in zip(frames['test']['proteins'], frames['test']['prop_annotations']) if p in keep
                  for t in sorted(set(a) & go_set)]
        shutil.copyfile(Path(root) / 'external/deepgose/data/go.obo', folder / 'go.obo')
    elif method == 'transfew':  # the ground-truth files exactly as given to CAFA-evaluator
        gt = transfew_groundtruth(root)
        labels = [tuple(line.rstrip('\n').split('\t')[:2]) for line in open(gt['Test_all'])]
        shutil.copyfile(gt['Test_novel'], folder / 'labels_test_novel.tsv')
        shutil.copyfile(tfew(root) / 'obo' / 'go-basic.obo', folder / 'go-basic.obo')
        shutil.copyfile(tfew(root) / 'test' / 'output_t1_t2' / 'ia_prop.txt', folder / 'ia_prop.txt')
    else:  # mlb.transform of mf_test2_go.txt; the evaluation adds is_a and part_of ancestors with the GO file
        truth = dpfunc_truth(root)
        labels = [(r['record'], t) for r in tests for t in sorted(truth[r['record']])]
        shutil.copyfile(Path(root) / 'datasets/ontology/go-basic.obo', folder / 'go.obo')
    write_tsv(folder / 'labels.tsv', ['record', 'go_id'], labels)
    log(f'testset: {len(tests):,} records, {len(labels):,} label rows, {len(vocabulary):,} terms')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('method', choices=list(PAPER))
    parser.add_argument('--testset', action='store_true', help='also write replicate/testset/')
    parser.add_argument('--paper-aupr', action='store_true',
                        help="also the paper's own protein-centric trapezoid AUPR (the papers' protocol check)")
    add_argument(parser)
    a = parser.parse_args()
    global PAPER_AUPR
    PAPER_AUPR = a.paper_aupr
    root = data_root(a.data)
    out = root / 'replicate' / 'results'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'paper.json').write_text(json.dumps(PAPER[a.method], indent=1) + '\n')
    globals()[a.method](root, out)
    if a.testset:
        testset(root, a.method)


if __name__ == '__main__':
    main()
