"""Convert the released ArchipelaGO dataset (Parquet, as downloaded into <data root>/release/) into dataset folders.

    python -m retrain.datasets.release [--data DIR]

Reads, under <data root>/release/:
  ontology/go-basic.obo.gz                     the GO release of every label
  go_terms/go_terms.parquet                    the 8,739 output terms; `column` is the output column
  training/train.parquet                       ArchipelaGO: protein_id, accessions, sequence, part (hq | selected),
                                               selection_rank, go_terms (propagated labels)
  training_annotations/train.parquet           asserted annotations: protein_id, go_id, evidence_class, evidence, source,
                                               assigned_by
  benchmark/{development,test}.parquet         protein_id, accessions, sequence, cluster, go_terms
  benchmark_annotations/{development,test}.parquet

Writes datasets/ontology (go-basic.obo, vocabulary.json), datasets/hq (the HQ rows of ArchipelaGO, which come first),
datasets/archipelago (columns part, selection_rank), datasets/development and datasets/test (column cluster). In
annotations.tsv the source of a GOA annotation is `GOA:<assigned_by>` (`GOA` for the older IBA/IEA panel labels without
an assigning group); every other source is written as released.

Check: for every dataset, propagating its annotations on the ontology (retrain.common.Ontology.label_matrix) must give
exactly the released labels; the conversion stops otherwise.
"""
import argparse
import gzip
import json
import shutil

import numpy as np
import pyarrow.parquet as pq
import scipy.sparse as sp

from retrain.common import Ontology, add_argument, data_root, log, move, write_dataset


def source_name(row):
    if row['source'] == 'GOA' and row['assigned_by']:
        return f"GOA:{row['assigned_by']}"
    return row['source']


def labels_from_terms(term_lists, ontology):
    rows, columns = [], []
    for i, terms in enumerate(term_lists):
        for t in terms:
            rows.append(i)
            columns.append(ontology.index[t])
    return sp.csr_matrix((np.ones(len(rows), np.uint8), (rows, columns)), shape=(len(term_lists), len(ontology.vocabulary)))


def annotations(table, wanted):
    """Annotation rows (seq_id, go_id, evidence, source) of the proteins in `wanted`."""
    out = []
    for r in table.to_pylist():
        if r['protein_id'] not in wanted:
            raise ValueError(f"annotation of a protein outside its table: {r['protein_id']}")
        out.append((r['protein_id'], r['go_id'], r['evidence'], source_name(r)))
    return out


def checked_labels(name, ids, term_lists, rows, ontology):
    """The released labels, after checking that the annotations propagate to exactly them."""
    released = labels_from_terms(term_lists, ontology)
    rebuilt, dropped = ontology.label_matrix([(h, g) for h, g, _, _ in rows], {h: i for i, h in enumerate(ids)})
    differ = (released.astype(np.int8) - rebuilt.astype(np.int8)).count_nonzero()
    if differ:
        raise ValueError(f'{name}: the annotations propagate to labels that differ from the released ones in {differ} pairs')
    log(f'{name}: annotations reproduce the {released.nnz:,} released labels ({len(dropped)} asserted pairs on terms that '
        'map to no live term)')
    return released


def ontology_folder(root, release):
    folder = root / 'datasets' / 'ontology'
    if folder.exists():
        raise FileExistsError(f'{folder} exists')
    terms = pq.read_table(release / 'go_terms' / 'go_terms.parquet', columns=['go_id', 'column']).to_pydict()
    vocabulary = [t for _, t in sorted(zip(terms['column'], terms['go_id']))]
    if sorted(terms['column']) != list(range(len(vocabulary))):
        raise ValueError('go_terms.parquet: columns are not 0 .. n-1')
    folder.mkdir(parents=True)
    with gzip.open(release / 'ontology' / 'go-basic.obo.gz', 'rb') as source, open(folder / 'go-basic.obo', 'wb') as target:
        shutil.copyfileobj(source, target)
    (folder / 'vocabulary.json').write_text(json.dumps(vocabulary) + '\n', encoding='utf-8')
    log(f'wrote datasets/ontology: {len(vocabulary):,} terms')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argument(parser)
    final = data_root(parser.parse_args().data)
    release = final / 'release'
    names = ('ontology', 'hq', 'archipelago', 'development', 'test')
    present = [n for n in names if (final / 'datasets' / n).exists()]
    if present:
        raise FileExistsError(f'datasets {present} exist under {final / "datasets"}')
    root = final / '.release.building'  # everything is built here, then moved into place: a failed run leaves nothing
    shutil.rmtree(root, ignore_errors=True)
    ontology_folder(root, release)
    ontology = Ontology(root)
    build(root, release, ontology)
    (final / 'datasets').mkdir(exist_ok=True)
    for name in names:
        move(root / 'datasets' / name, final / 'datasets' / name)
    shutil.rmtree(root)
    log(f'datasets {", ".join(names)} written under {final / "datasets"}')


def build(root, release, ontology):
    """The four datasets of the release, under <root>/datasets/."""
    training = pq.read_table(release / 'training' / 'train.parquet').to_pydict()
    ids = training['protein_id']
    parts = training['part']
    n_hq = parts.count('hq')
    if parts != ['hq'] * n_hq + ['selected'] * (len(ids) - n_hq):
        raise ValueError('training table: rows are not the HQ proteins followed by the selected proteins')
    ranks = training['selection_rank']
    if ranks[n_hq:] != list(range(1, len(ids) - n_hq + 1)):
        raise ValueError('training table: selected proteins are not in ranking order')
    rows = annotations(pq.read_table(release / 'training_annotations' / 'train.parquet'), set(ids))
    labels = checked_labels('archipelago', ids, training['go_terms'], rows, ontology)
    hq_ids = set(ids[:n_hq])
    write_dataset(root, 'hq', ids[:n_hq], training['accessions'][:n_hq], training['sequence'][:n_hq],
                  [r for r in rows if r[0] in hq_ids], labels[:n_hq], ontology)
    write_dataset(root, 'archipelago', ids, training['accessions'], training['sequence'], rows, labels, ontology,
                  columns=dict(part=parts, selection_rank=['' if r is None else r for r in ranks]))

    for split in ('development', 'test'):
        table = pq.read_table(release / 'benchmark' / f'{split}.parquet').to_pydict()
        rows = annotations(pq.read_table(release / 'benchmark_annotations' / f'{split}.parquet'), set(table['protein_id']))
        labels = checked_labels(split, table['protein_id'], table['go_terms'], rows, ontology)
        write_dataset(root, split, table['protein_id'], table['accessions'], table['sequence'], rows, labels, ontology,
                      columns=dict(cluster=table['cluster']))


if __name__ == '__main__':
    main()
