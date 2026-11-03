"""A published method's own split as a data root of this repository's code.

    python -m retrain.replicate.build {deepgose,transfew,dpfunc} --data DIR

DIR is a data root for this method alone; external/ holds the authors' release (README, External inputs):
  DeepGO-SE  external/deepgose/data/ (training-data.tar.gz unpacked): mf/{train,valid,test}_data.pkl, 38,533 / 1,901 /
             2,845 records; output terms mf/terms.pkl (6,851); labels = each record's prop_annotations within the terms
             (the official load_data); GO = the archive's go.obo (releases/2021-11-16); axioms = its go-plus.norm.
  TransFew   external/transfew/TFewData (the extracted release): train 70,773 and validation 7,864 records (the released
             label vectors over mf/sorted_terms, 7,224 terms, already propagated; retrain.transfew.original.records), test
             = Test_all, the 1,112 Swiss-Prot proteins of test/output_t1_t2/groundtruths/full/mf_swissprot.tsv (Test_novel:
             the 858 of groundtruths/seq_ID_30); GO = TFewData/obo/go-basic.obo (releases/2023-01-01).
  DPFunc     external/dpfunc/dpfunc_published_dataset.zip: the used records of mf_{train,test1,test2}, 31,463 / 682 / 401
             (retrain.dpfunc.original.records); output terms = the 6,086 terms of the training labels (what DPFunc's get_mlb
             fits; checked against the repository's mlb/mf_go.mlb in external/dpfunc/DPFunc, commit c5faf4df); labels =
             the supplied (already propagated) terms within them, as mlb.transform; GO = external/go/2022-12-04/go.obo
             (http://release.geneontology.org/2022-12-04/ontology/go.obo; the release does not include its GO file, and
             2022-12-04 is the only release from May 2020 to April 2023 under which all its labels are closed).
Records with the same sequence are one protein within a split (labels united); identical sequences in different splits stay
in each (the authors' split as released; train with --allow-shared-proteins). Vocabulary terms that are not live in the GO
file are listed in datasets/ontology/not_live.json (DPFunc: 77). Test labels written here are for reference: score.py
rebuilds the truth with each paper's own code.

Writes datasets/{ontology, original_<method>, development, test}, replicate/records.tsv (split, record, accession, seq_id[,
novel]) and, for DeepGO-SE, inputs/deepgose/axioms/go-plus.norm. Then compute the method's inputs with its own stages
(README) on these three datasets.
"""
import argparse
import collections
import json
import shutil
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from retrain.common import Ontology, add_argument, data_root, log, seq_id, write_dataset

DATASET = {'deepgose': 'original_deepgose', 'transfew': 'original_transfew', 'dpfunc': 'original_dpfunc'}
SPLIT_NAME = dict(train=None, valid='development', test='test')


def external(root, *parts):
    path = Path(root, 'external', *parts)
    if not path.exists():
        raise SystemExit(f'{path} is missing (README, External inputs)')
    return path


def ontology(root, obo, vocabulary):
    """datasets/ontology: the authors' GO file (as go-basic.obo), their vocabulary, and not_live.json if needed."""
    from archipelago.go import read_obo
    folder = Path(root) / 'datasets' / 'ontology'
    if folder.exists():
        raise FileExistsError(f'{folder} exists; build a replication root once')
    folder.mkdir(parents=True)
    shutil.copyfile(obo, folder / 'go-basic.obo')
    (folder / 'vocabulary.json').write_text(json.dumps(list(vocabulary)) + '\n')
    terms = read_obo(folder / 'go-basic.obo')
    not_live = [t for t in vocabulary if t not in terms or terms[t]['obsolete']]
    if not_live:
        (folder / 'not_live.json').write_text(json.dumps(not_live) + '\n')
    log(f'ontology: {obo}, {len(vocabulary):,} terms, {len(not_live)} not live in it')
    return Ontology(root)


def write(root, method, rows, onto, extra=()):
    """rows: dict(split, record, accession, sequence, terms[, extra columns]). Writes the three datasets (one protein per
    distinct sequence within a split) and replicate/records.tsv."""
    for split, name in (('train', DATASET[method]), ('valid', 'development'), ('test', 'test')):
        proteins = {}
        for r in rows:
            if r['split'] != split:
                continue
            p = proteins.setdefault(seq_id(r['sequence']), dict(sequence=r['sequence'].upper(), terms=set(), records=[],
                                                                  accessions=[]))
            p['terms'] |= set(r['terms'])
            p['records'].append(r['record'])
            p['accessions'] += [a for a in r['accession'].split(',') if a and a not in p['accessions']]
        ids = list(proteins)
        labels = sp.lil_matrix((len(ids), len(onto.vocabulary)), dtype=np.uint8)
        for i, h in enumerate(ids):
            for t in proteins[h]['terms']:
                labels[i, onto.index[t]] = 1
        annotations = [(h, t, 'released', f'{method} {split}') for h in ids for t in sorted(proteins[h]['terms'])]
        write_dataset(root, name, ids, [proteins[h]['accessions'] for h in ids], [proteins[h]['sequence'] for h in ids],
                      annotations, labels.tocsr(), onto, columns=dict(records=[','.join(proteins[h]['records']) for h in ids]))
    folder = Path(root) / 'replicate'
    folder.mkdir(exist_ok=True)
    with open(folder / 'records.tsv', 'w') as handle:
        handle.write('\t'.join(['split', 'record', 'accession', 'seq_id', *extra]) + '\n')
        handle.writelines('\t'.join([r['split'], r['record'], r['accession'], seq_id(r['sequence']),
                                     *[str(r.get(c, '')) for c in extra]]) + '\n' for r in rows)
    log(f"{method}: records {dict(collections.Counter(r['split'] for r in rows))}")


def deepgose(root):
    import pandas as pd
    data = external(root, 'deepgose', 'data')
    terms = list(pd.read_pickle(data / 'mf' / 'terms.pkl')['gos'].values.flatten())
    vocabulary, rows = set(terms), []
    for split, expected in (('train', 38533), ('valid', 1901), ('test', 2845)):
        df = pd.read_pickle(data / 'mf' / f'{split}_data.pkl')
        if len(df) != expected:
            raise ValueError(f'{split}_data.pkl: {len(df)} records, expected {expected}')
        rows += [dict(split=split, record=r.proteins, sequence=r.sequences, terms=set(r.prop_annotations) & vocabulary,
                      accession=','.join(a.strip() for a in r.accessions.split(';') if a.strip())) for r in df.itertuples()]
    write(root, 'deepgose', rows, ontology(root, data / 'go.obo', terms))
    axioms = Path(root) / 'inputs' / 'deepgose' / 'axioms'
    axioms.mkdir(parents=True)
    shutil.copyfile(data / 'go-plus.norm', axioms / 'go-plus.norm')


def transfew(root):
    from retrain.transfew.original import index_pickle, records
    release = external(root, 'transfew', 'TFewData')
    terms = index_pickle(release, 'sorted_terms')
    novel = {line.split('\t')[0] for line in open(release / 'test/output_t1_t2/groundtruths/seq_ID_30/mf_swissprot.tsv')}
    split_of = dict(train='train', validation='valid', test_swissprot='test')  # test_trembl: not in the paper's table
    rows = [dict(split=split_of[r['split']], record=r['record'], accession=r['accession'], sequence=r['sequence'],
                 terms=set(r['terms']) & set(terms), novel=int(r['accession'] in novel) if r['split'] == 'test_swissprot' else '')
            for r in records(release) if r['split'] in split_of]
    counts = collections.Counter(r['split'] for r in rows)
    if counts != dict(train=70773, valid=7864, test=1112):
        raise ValueError(f'TransFew records {dict(counts)}')
    write(root, 'transfew', rows, ontology(root, release / 'obo' / 'go-basic.obo', terms), extra=('novel',))


def dpfunc(root):
    import joblib
    from retrain.dpfunc.original import records
    recs = records(root)
    classes = sorted(set().union(*[t for _, s, _, t in recs if s == 'train']))
    released = list(joblib.load(external(root, 'dpfunc', 'DPFunc', 'mlb', 'mf_go.mlb')).classes_)
    if released != classes:
        raise ValueError(f'mlb/mf_go.mlb ({len(released)} classes) differs from the training terms ({len(classes)})')
    rows = [dict(split=split, record=name, accession='', sequence=sequence, terms=set(terms) & set(classes))
            for name, split, sequence, terms in recs]
    write(root, 'dpfunc', rows, ontology(root, external(root, 'go', '2022-12-04', 'go.obo'), classes))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('method', choices=list(DATASET))
    add_argument(parser)
    a = parser.parse_args()
    globals()[a.method](data_root(a.data))


if __name__ == '__main__':
    main()
