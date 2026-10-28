"""DeepGO-SE's original training data (datasets/original_deepgose), rebuilt from the authors' processed archive.

Input: the deepgo2 training data, https://bio2vec.net/data/deepgo/deepgo2/training-data.tar.gz (SHA-256
80cc5446ef2d32cfccd19ff6a634cb6c2946370caf0ec7b9fd1bfb3c2dc07a44), unpacked so that
<root>/external/deepgose/data/mf/{train,valid,test}_data.pkl exist (43,279 records, the MF count of the paper's Table 6).
Each record is a Swiss-Prot entry with `annotations` (GO id|evidence code) and `exp_annotations` (the GO ids with an
experimental code: EXP IDA IPI IMP IGI IEP TAS IC HTP HDA HMP HGI HEP, deepgo/utils.py EXP_CODES).

Rules:
- The three published splits are all training data (the authors' validation and test proteins are not needed here).
- One protein per distinct sequence; records with the same sequence are merged (union of annotations).
- Screen: proteins identical to, or with a DIAMOND hit to, a development or test protein are removed (retrain.common.screen).
- Labels: each record's experimental annotations (not the supplied `prop_annotations`, which were propagated on the
  authors' 2021 ontology over every relationship type), mapped to live terms of datasets/ontology (alt_id -> primary id;
  an obsolete term with a single live replaced_by -> that term; other obsolete terms dropped), MF terms kept, propagated
  over is_a and restricted to the vocabulary.
- Inclusion: the authors' rule (an MF term among the experimental annotations), applied on this ontology: proteins
  without any live MF experimental annotation are left out.

    python -m retrain.deepgose.original [--threads 16]

Writes datasets/original_deepgose and, in inputs/deepgose/original/, removed.tsv (proteins removed by the screen) and
without_live_mf.tsv (proteins left out by the inclusion rule). CPU only; minutes.
"""
import argparse
import collections
from pathlib import Path

import pandas as pd

from retrain.common import Ontology, add_argument, data_root, log, screen, seq_id, write_dataset

EXP_CODES = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'}
SPLITS = ['train', 'valid', 'test']


def records(source):
    """Every published MF record, in train, valid, test order, with its experimental GO ids and their codes."""
    out = []
    for split in SPLITS:
        for row in pd.read_pickle(Path(source) / f'{split}_data.pkl').itertuples():
            codes = collections.defaultdict(set)
            for item in row.annotations:
                go_id, code = item.split('|')
                if code in EXP_CODES:
                    codes[go_id].add(code)
            if set(codes) != set(row.exp_annotations):
                raise ValueError(f'{split}:{row.proteins}: exp_annotations differ from the experimental annotations')
            if row.sequences != row.sequences.upper():
                raise ValueError(f'{split}:{row.proteins}: lower-case sequence')
            out.append(dict(name=f'{split}:{row.proteins}', sequence=row.sequences, codes=codes,
                            accessions=[a.strip() for a in row.accessions.split(';') if a.strip()]))
    return out


def merge(recs):
    """{seq_id: protein} with one protein per distinct sequence, in order of first appearance."""
    proteins = {}
    for r in recs:
        p = proteins.setdefault(seq_id(r['sequence']), dict(sequence=r['sequence'], records=[], accessions=[],
                                                             codes=collections.defaultdict(set)))
        p['records'].append(r['name'])
        p['accessions'] += [a for a in r['accessions'] if a not in p['accessions']]
        for go_id, codes in r['codes'].items():
            p['codes'][go_id] |= codes
    return proteins


def mf_assertions(ontology, codes):
    """[(asserted id, live MF term, evidence codes)] of a protein's experimental annotations."""
    out = []
    for go_id, evidence in sorted(codes.items()):
        term = ontology.map_term(go_id)
        if term is not None and ontology.terms[term]['namespace'] == 'molecular_function':
            out.append((go_id, term, evidence))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--threads', type=int, default=8, help='DIAMOND threads')
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    work = root / 'inputs' / 'deepgose' / 'original'
    ontology = Ontology(root)
    recs = records(root / 'external' / 'deepgose' / 'data' / 'mf')
    proteins = merge(recs)
    log(f'{len(recs):,} records, {len(proteins):,} distinct sequences')

    exact, hits = screen(root, {h: p['sequence'] for h, p in proteins.items()}, work / 'screen', a.threads)
    removed = exact | hits
    screened = [h for h in proteins if h not in removed]
    mf = {h: mf_assertions(ontology, proteins[h]['codes']) for h in screened}
    kept = [h for h in screened if mf[h]]
    left_out = [h for h in screened if not mf[h]]
    log(f'{len(removed):,} removed by the screen, {len(left_out):,} without a live MF experimental annotation, '
        f'{len(kept):,} kept')

    row_of = {h: i for i, h in enumerate(kept)}
    pairs, rows = [], []
    for h in kept:
        source = 'DeepGO-SE MF ' + ';'.join(proteins[h]['records'])
        for go_id, term, evidence in mf[h]:
            pairs.append((h, term))
            rows += [(h, term, code, source + ('' if term == go_id else f' (mapped from {go_id})'))
                     for code in sorted(evidence)]
    labels, dropped = ontology.label_matrix(pairs, row_of)
    if dropped or (labels.getnnz(axis=1) == 0).any():
        raise ValueError('a kept protein has an unmapped term or no label in the vocabulary')
    write_dataset(root, 'original_deepgose', kept, [proteins[h]['accessions'] for h in kept],
                  [proteins[h]['sequence'] for h in kept], rows, labels, ontology,
                  columns={'records': [','.join(proteins[h]['records']) for h in kept]})
    with open(work / 'removed.tsv', 'w') as f:
        f.write('seq_id\trecords\treason\n')
        f.writelines(f"{h}\t{','.join(proteins[h]['records'])}\t{'identical' if h in exact else 'diamond'}\n"
                     for h in proteins if h in removed)
    with open(work / 'without_live_mf.tsv', 'w') as f:
        f.write('seq_id\trecords\texperimental_annotations\n')
        f.writelines(f"{h}\t{','.join(proteins[h]['records'])}\t{','.join(sorted(proteins[h]['codes']))}\n"
                     for h in left_out)


if __name__ == '__main__':
    main()
