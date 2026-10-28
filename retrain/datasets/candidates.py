"""datasets/candidates: the computational candidates of the coverage ranking, from their release files.

    python -m retrain.datasets.candidates [--data DIR]

The candidates are the 8,952,628 TrEMBL proteins with computational molecular-function annotations only that remain
after one random representative is taken per 50%-identity cluster (DIAMOND DeepClust), relatives of the benchmark are
removed (exact sequence or DIAMOND hit, E <= 1e-3, either direction) and proteins with a label outside the vocabulary are
left out. Release files (to be released with the dataset), under <data root>/release/candidates/:

  proteins.parquet     protein_id (SHA-256 of the sequence), accessions (list of UniProt accessions), sequence,
                       go_terms (list: propagated labels, vocabulary terms only). Rows in the released pool order, the
                       order that the paper's random control orders permute (retrain.datasets.selected); about 3 GB.
  annotations.parquet  protein_id, go_id, evidence, assigned_by: the computational GOA annotations behind the labels;
                       about 20 million rows, 0.3 GB.

Writes datasets/candidates (proteins.tsv, sequences.fasta, annotations.tsv, labels.npz) in the released row order, after
checking that every sequence hashes to its protein_id, that every candidate has a label, and that the annotations
propagate to exactly the released labels. About 15 minutes and 16 GB of memory.
"""
import argparse
import csv
import shutil

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import scipy.sparse as sp

from retrain.common import Ontology, add_argument, data_root, log, move, seq_id


def term_matrix(table, ontology):
    """(released labels CSR, protein ids) from the proteins table's go_terms lists."""
    ids = table.column('protein_id').combine_chunks()
    terms = table.column('go_terms').combine_chunks()
    lengths = pc.list_value_length(terms).to_numpy(zero_copy_only=False)
    columns = pc.index_in(pc.list_flatten(terms), value_set=pa.array(ontology.vocabulary))
    if columns.null_count:
        raise ValueError('a candidate label is not a vocabulary term')
    if (lengths == 0).any():
        raise ValueError(f'{int((lengths == 0).sum())} candidates have no label')
    indptr = np.concatenate([[0], np.cumsum(lengths)])
    labels = sp.csr_matrix((np.ones(indptr[-1], np.uint8), columns.to_numpy(), indptr),
                           shape=(len(ids), len(ontology.vocabulary)))
    labels.sum_duplicates()
    labels.data[:] = 1
    return labels, ids


def propagated(annotations, ids, ontology):
    """Labels from the annotations: (proteins x asserted terms) @ (asserted terms x their vocabulary ancestors) > 0."""
    rows = pc.index_in(annotations.column('protein_id'), value_set=ids)
    if rows.null_count:
        raise ValueError('an annotation names a protein that is not a candidate')
    asserted = annotations.column('go_id').combine_chunks().dictionary_encode()
    terms = asserted.dictionary.to_pylist()
    closure = sp.lil_matrix((len(terms), len(ontology.vocabulary)), dtype=np.uint8)
    unmapped = 0
    for k, go_id in enumerate(terms):
        live = ontology.map_term(go_id)
        if live is None:
            unmapped += 1
            continue
        for ancestor in ontology.ancestors(live):
            if ancestor in ontology.index:
                closure[k, ontology.index[ancestor]] = 1
    a = sp.csr_matrix((np.ones(len(rows), np.float32), (rows.to_numpy(zero_copy_only=False), asserted.indices.to_numpy())),
                      shape=(len(ids), len(terms)))
    out = (a @ closure.tocsr().astype(np.float32)).tocsr()
    out.data[:] = 1
    return out.astype(np.uint8), unmapped


def checked(source, ontology):
    """(labels CSR, protein ids, annotations table) of a release folder holding proteins.parquet and annotations.parquet
    in this format, after checking that the annotations propagate to exactly the released labels."""
    labels, ids = term_matrix(pq.read_table(source / 'proteins.parquet', columns=['protein_id', 'go_terms']), ontology)
    log(f'{source.name}: {len(ids):,} proteins, {labels.nnz:,} labels')
    annotations = pq.read_table(source / 'annotations.parquet')
    rebuilt, unmapped = propagated(annotations, ids, ontology)
    differ = (rebuilt.astype(np.int8) - labels.astype(np.int8)).count_nonzero()
    if differ:
        raise ValueError(f'the annotations propagate to labels that differ from the released ones in {differ} pairs')
    log(f'annotations ({annotations.num_rows:,} rows, {unmapped} terms without a live id) reproduce the labels')
    return labels, ids, annotations


def annotation_row(h, go_id, evidence, assigned_by):
    """An annotations.tsv row (seq_id, go_id, evidence, source) of a released computational GOA annotation."""
    return h, go_id, evidence or '', f'GOA:{assigned_by}' if assigned_by else 'GOA'


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argument(parser)
    root = data_root(parser.parse_args().data)
    source, folder = root / 'release' / 'candidates', root / 'datasets' / 'candidates'
    if folder.exists():
        raise FileExistsError(f'{folder} exists')
    labels, ids, annotations = checked(source, Ontology(root))

    building = folder.with_name('.candidates.building')
    shutil.rmtree(building, ignore_errors=True)
    building.mkdir(parents=True)
    with open(building / 'proteins.tsv', 'w', newline='', encoding='utf-8') as table, \
            open(building / 'sequences.fasta', 'w', encoding='utf-8') as fasta:
        writer = csv.writer(table, delimiter='\t', lineterminator='\n')
        writer.writerow(['seq_id', 'accessions'])
        done = 0
        for batch in pq.ParquetFile(source / 'proteins.parquet').iter_batches(
                batch_size=500_000, columns=['protein_id', 'accessions', 'sequence']):
            for h, accessions, sequence in zip(*(batch.column(k).to_pylist() for k in ('protein_id', 'accessions',
                                                                                       'sequence'))):
                if h != ids[done].as_py() or seq_id(sequence) != h:
                    raise ValueError(f'row {done}: {h} out of order or not the SHA-256 of its sequence')
                writer.writerow([h, ','.join(accessions or [])])
                fasta.write(f'>{h}\n{sequence.upper()}\n')
                done += 1
            log(f'proteins written: {done:,} of {len(ids):,}')
    with open(building / 'annotations.tsv', 'w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['seq_id', 'go_id', 'evidence', 'source'])
        for batch in annotations.to_batches(1_000_000):
            columns = [batch.column(k).to_pylist() for k in ('protein_id', 'go_id', 'evidence', 'assigned_by')]
            writer.writerows(annotation_row(*row) for row in zip(*columns))
    sp.save_npz(building / 'labels.npz', labels, compressed=True)
    move(building, folder)
    log(f'wrote {folder}')


if __name__ == '__main__':
    main()
