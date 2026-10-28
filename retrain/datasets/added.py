"""Training sets made of HQ followed by added computational proteins: reading the inputs they come from and writing them.

Every such set has the HQ rows first (datasets/hq: proteins, labels and annotations unchanged), then the added proteins in
the order of their selection rule. An added protein's labels are its computational GO annotations propagated over is_a
within the vocabulary.
"""
import csv

import numpy as np
import pyarrow as pa
import pyarrow.csv as pacsv
import scipy.sparse as sp

from retrain.common import load, log, read_annotations, write_dataset

RANKINGS = {  # the paper's coverage rankings (release files) and their lengths
    'coverage_ranking.tsv': 600_000,           # with the rare-function filter: ArchipelaGO's 150,000 and the dose sets
    'coverage_ranking_unfiltered.tsv': 500_000,  # without the filter: the coverage_unfiltered control (first 165,000)
}


def column(path, name):
    """One column of a large tab-separated table (all values as strings), read with pyarrow."""
    table = pacsv.read_csv(path, parse_options=pacsv.ParseOptions(delimiter='\t'),
                           convert_options=pacsv.ConvertOptions(column_types={name: pa.string()}, include_columns=[name]))
    return table.column(name).to_pylist()


def sequences(path, wanted):
    """{id: sequence} of the wanted ids of a FASTA file, streamed (the candidates' file has 9 million sequences)."""
    out, keep, name, parts = {}, False, None, []
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('>'):
                if keep:
                    out[name] = ''.join(parts)
                name = line[1:].split()[0]
                keep, parts = name in wanted, []
            elif keep:
                parts.append(line.strip())
    if keep:
        out[name] = ''.join(parts)
    missing = wanted - set(out)
    if missing:
        raise KeyError(f'{path}: {len(missing)} wanted proteins missing, e.g. {next(iter(missing))}')
    return out


def rows_of(path, wanted):
    """Rows (seq_id, go_id, evidence, source) of an annotations.tsv whose protein is wanted. Streamed rather than read
    with retrain.common.read_annotations: the candidates' file has tens of millions of rows."""
    out = []
    with open(path, newline='', encoding='utf-8') as handle:
        for r in csv.DictReader(handle, delimiter='\t'):
            if r['seq_id'] in wanted:
                out.append((r['seq_id'], r['go_id'], r['evidence'], r['source']))
    return out


class Candidates:
    """datasets/candidates: the computational candidates that the coverage ranking and the control sets choose from (one
    representative per 50%-identity cluster of computationally annotated TrEMBL proteins, related proteins of the benchmark
    removed, labels inside the vocabulary; retrain.datasets.candidates). Its row order is the released pool order, which
    the paper's random orders permute."""

    def __init__(self, root):
        self.folder = root / 'datasets' / 'candidates'
        self.ids = column(self.folder / 'proteins.tsv', 'seq_id')
        self.row = {h: i for i, h in enumerate(self.ids)}
        if len(self.row) != len(self.ids):
            raise ValueError(f'{self.folder}: repeated seq_id')
        self.labels = sp.load_npz(self.folder / 'labels.npz').tocsr()
        if self.labels.shape[0] != len(self.ids):
            raise ValueError(f'{self.folder}: labels have {self.labels.shape[0]} rows for {len(self.ids)} proteins')
        log(f'candidates: {len(self.ids):,} proteins')

    def rows(self, ids):
        missing = [h for h in ids if h not in self.row]
        if missing:
            raise KeyError(f'{len(missing)} proteins are not candidates, e.g. {missing[0]}')
        return [self.row[h] for h in ids]

    def details(self, ids):
        """(accessions, sequences, annotation rows) of these candidates."""
        wanted = set(ids)
        accessions = dict(zip(self.ids, column(self.folder / 'proteins.tsv', 'accessions')))
        found = sequences(self.folder / 'sequences.fasta', wanted)
        return ([accessions[h].split(',') if accessions[h] else [] for h in ids], [found[h] for h in ids],
                rows_of(self.folder / 'annotations.tsv', wanted))


def rare_function(labels, hq_counts, below=1000):
    """The rare-function filter: True for a candidate with at least one label held by fewer than `below` HQ proteins."""
    labels = sp.csr_matrix(labels)
    rarest = np.full(labels.shape[0], np.iinfo(np.int64).max)
    has = np.diff(labels.indptr) > 0
    rarest[has] = np.minimum.reduceat(hq_counts[labels.indices], labels.indptr[:-1][has])
    return rarest < below


def ranking(root, name):
    """seq_ids of a coverage ranking (rank, seq_id, score), in rank order. The paper's ranking, released as
    release/selection/<name>, is used when present (its length checked); otherwise selection/<name> from
    scripts/rank_coverage.py, with a warning: a recomputed ranking can order exactly tied candidates differently, and
    only ArchipelaGO's first 150,000 picks are checked against the release."""
    released = root / 'release' / 'selection' / name
    path = released if released.is_file() else root / 'selection' / name
    with open(path, newline='', encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if [int(r['rank']) for r in rows] != list(range(1, len(rows) + 1)):
        raise ValueError(f'{path}: ranks are not 1, 2, 3, ...')
    if path == released and len(rows) != RANKINGS[name]:
        raise ValueError(f'{path}: {len(rows):,} rows; the paper\'s ranking has {RANKINGS[name]:,}')
    if path != released:
        log(f'WARNING: {released} not found; using the recomputed {path}. Sets beyond ArchipelaGO built from it can differ '
            "from the paper's (tie order, and embeddings of non-standard residues).")
    return [r['seq_id'] for r in rows]


def write_added(root, name, ontology, ids, ranks, accessions, seqs, rows, labels, part, rank_column):
    """Write <name> = datasets/hq followed by the added proteins. Columns: part = hq | <part>, and <rank_column> = an added
    protein's 1-based position in its selection order (proteins skipped by a filter keep their positions), empty for HQ."""
    hq = load(root, 'hq')
    held = set(load(root, 'development', labels=False).seq_ids) | set(load(root, 'test', labels=False).seq_ids)
    overlap = set(ids) & (set(hq.seq_ids) | held)
    if overlap:
        raise ValueError(f'{name}: {len(overlap)} added proteins are HQ, development or test proteins')
    labels = sp.csr_matrix(labels)
    if (np.diff(labels.indptr) == 0).any():
        raise ValueError(f'{name}: an added protein has no label')
    n_hq = len(hq.seq_ids)
    write_dataset(root, name, hq.seq_ids + list(ids),
                  [a.split(',') if a else [] for a in hq.accessions] + list(accessions), hq.sequences + list(seqs),
                  read_annotations(hq.folder) + list(rows), sp.vstack([hq.labels, labels]).tocsr(), ontology,
                  columns={'part': ['hq'] * n_hq + [part] * len(ids), rank_column: [''] * n_hq + [str(k) for k in ranks]})
    log(f'{name}: {n_hq:,} HQ + {len(ids):,} added proteins with {labels.nnz:,} labels')
