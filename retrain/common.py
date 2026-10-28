"""Shared by every stage of the retraining pipeline: dataset folders, the GO ontology, fit folders and the homology screen.

Everything lives under the data root (retrain.paths; layout in data/README.md):

    datasets/<name>/    proteins.tsv, sequences.fasta, annotations.tsv, labels.npz   (ontology/: go-basic.obo, vocabulary.json)
    inputs/<kind>/...   derived inputs (embeddings, features, structures, ...)
    runs/<method>/<training set>/seed<k>/   dev_scores.npy, test_scores.npy, record.json

A protein's identifier (seq_id) is the SHA-256 hex digest of its upper-case sequence. Labels are a CSR uint8 matrix
(proteins x vocabulary terms, scipy.sparse.save_npz), propagated over is_a within the vocabulary.
"""
import csv
import hashlib
import json
import os
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from archipelago import io
from archipelago.go import Propagation, read_obo
from retrain.paths import add_argument, data_root  # noqa: F401  (re-exported for the stages)

SPLITS = {'dev': 'development', 'test': 'test'}  # file prefix of a fit's scores -> dataset name
DIAMOND_DBSIZE = 95_692_462  # fixed database size of every homology screen in the paper (E-values comparable across runs)


def log(*args):
    print(time.strftime('%Y-%m-%d %H:%M:%S'), *args, flush=True)


def seq_id(sequence):
    return hashlib.sha256(sequence.upper().encode()).hexdigest()


def read_fasta(path):
    """(ids, sequences) in file order; the id is the first word after '>'."""
    records = io.read_fasta(path)
    return list(records), list(records.values())


def write_fasta(path, items):
    with open(path, 'w', encoding='utf-8') as handle:
        handle.writelines(f'>{name}\n{sequence}\n' for name, sequence in items)


# ---------------------------------------------------------------------------------------------------------- datasets
@dataclass
class Dataset:
    name: str
    folder: Path
    seq_ids: list
    accessions: list
    sequences: list
    labels: object = None  # scipy CSR uint8 (proteins x terms), or None
    columns: dict = field(default_factory=dict)  # further proteins.tsv columns, as lists


def load(root, name, labels=True, scoring=False, sequences=True):
    """The dataset `name` under <root>/datasets/. The test labels are read only for scoring (scoring=True), so no
    training or selection code can use them by accident. sequences=False skips sequences.fasta (large pools)."""
    if name == 'test' and labels and not scoring:
        raise PermissionError('test labels are read only for scoring: load(..., scoring=True)')
    folder = Path(root) / 'datasets' / name
    rows = io.read_table(folder / 'proteins.tsv')
    ids = [r['seq_id'] for r in rows]
    seqs = None
    if sequences:
        fasta = io.read_fasta(folder / 'sequences.fasta')
        if list(fasta) != ids:
            raise ValueError(f'{folder}: sequences.fasta is not in proteins.tsv order')
        seqs = list(fasta.values())
    y = None
    if labels:
        y = sp.load_npz(folder / 'labels.npz').tocsr()
        if y.shape[0] != len(ids):
            raise ValueError(f'{folder}: labels have {y.shape[0]} rows for {len(ids)} proteins')
    columns = {k: [r[k] for r in rows] for k in (rows[0] if rows else {}) if k not in ('seq_id', 'accessions')}
    return Dataset(name, folder, ids, [r['accessions'] for r in rows], seqs, y, columns)


def write_dataset(root, name, seq_ids, accessions, sequences, annotations, labels, ontology, columns=None):
    """Write <root>/datasets/<name> (which must not exist). annotations: (seq_id, go_id, evidence, source) rows;
    accessions: lists; columns: {name: list} written after seq_id and accessions. Labels must be closed under is_a."""
    folder = Path(root) / 'datasets' / name
    if folder.exists():
        raise FileExistsError(f'{folder} exists; write a dataset once (remove it to rebuild)')
    n = len(seq_ids)
    if len(set(seq_ids)) != n or len(sequences) != n or len(accessions) != n:
        raise ValueError(f'{name}: identifiers repeat, or ids, sequences and accessions differ in number')
    wrong = [h for h, s in zip(seq_ids, sequences) if seq_id(s) != h]
    if wrong:
        raise ValueError(f'{name}: {len(wrong)} sequences do not hash to their seq_id, e.g. {wrong[0]}')
    labels = sp.csr_matrix(labels, dtype=np.uint8)
    labels.sort_indices()
    if labels.shape != (n, len(ontology.vocabulary)) or (labels.data != 1).any():
        raise ValueError(f'{name}: labels of shape {labels.shape} or with values other than 1')
    ontology.assert_closed(labels)
    rows = sorted(set(annotations))
    known = set(seq_ids)
    if any(r[0] not in known for r in rows):
        raise ValueError(f'{name}: an annotation names a protein outside the dataset')
    building = folder.with_name(f'.{name}.building')
    shutil.rmtree(building, ignore_errors=True)
    building.mkdir(parents=True)
    extra = list(columns or {})
    with open(building / 'proteins.tsv', 'w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['seq_id', 'accessions', *extra])
        for i, h in enumerate(seq_ids):
            writer.writerow([h, ','.join(accessions[i]), *[columns[c][i] for c in extra]])
    write_fasta(building / 'sequences.fasta', zip(seq_ids, sequences))
    with open(building / 'annotations.tsv', 'w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['seq_id', 'go_id', 'evidence', 'source'])
        writer.writerows(rows)
    sp.save_npz(building / 'labels.npz', labels, compressed=True)
    move(building, folder)
    log(f'wrote datasets/{name}: {n:,} proteins, {labels.nnz:,} propagated labels, {len(rows):,} annotations')
    return folder


def move(source, target, attempts=10):
    """Rename a just-written folder into place. On Windows a virus scanner can hold new files for a moment, which makes
    the rename fail with 'Access is denied'; it is retried for a few seconds before failing."""
    for attempt in range(attempts):
        try:
            return Path(source).rename(target)
        except PermissionError:
            if attempt == attempts - 1:
                raise
            time.sleep(0.5)


def read_annotations(folder):
    """The asserted annotations of a dataset folder: (seq_id, go_id, evidence, source) rows of its annotations.tsv."""
    with open(Path(folder) / 'annotations.tsv', newline='', encoding='utf-8') as handle:
        return [(r['seq_id'], r['go_id'], r['evidence'], r['source']) for r in csv.DictReader(handle, delimiter='\t')]


def label_counts(labels):
    labels = sp.csr_matrix(labels)
    return dict(proteins=labels.shape[0], propagated_pairs=int(labels.nnz),
                terms_with_positive=int((labels.getnnz(axis=0) > 0).sum()))


# ---------------------------------------------------------------------------------------------------------- ontology
class Ontology:
    """The GO release of every label (datasets/ontology/go-basic.obo) and the output vocabulary (vocabulary.json)."""

    def __init__(self, root):
        folder = Path(root) / 'datasets' / 'ontology'
        self.obo = folder / 'go-basic.obo'
        self.terms = read_obo(self.obo)
        self.vocabulary = json.loads((folder / 'vocabulary.json').read_text(encoding='utf-8'))
        self.index = {t: i for i, t in enumerate(self.vocabulary)}
        if len(self.index) != len(self.vocabulary):
            raise ValueError('vocabulary.json has duplicate terms')
        bad = [t for t in self.vocabulary if t not in self.terms or self.terms[t]['obsolete']]
        if bad:
            raise ValueError(f'{len(bad)} vocabulary terms are missing or obsolete in go-basic.obo, e.g. {bad[0]}')
        self.alt = {a: t for t, d in self.terms.items() for a in d['alt_id']}
        self.parents = [sorted({self.index[p] for p in self.terms[t]['is_a'] if p in self.index}) for t in self.vocabulary]
        self._propagation = Propagation(self.parents)
        self._ancestors = {}

    def map_term(self, go_id):
        """The live id of go_id: an alt_id maps to its primary id; an obsolete term with exactly one replaced_by that is
        live maps to it (GO's designated replacement; `consider` is not followed); None for anything else."""
        go_id = self.alt.get(go_id, go_id)
        term = self.terms.get(go_id)
        if term is None:
            return None
        if not term['obsolete']:
            return go_id
        if len(term['replaced_by']) == 1:
            replacement = self.alt.get(term['replaced_by'][0], term['replaced_by'][0])
            if replacement in self.terms and not self.terms[replacement]['obsolete']:
                return replacement
        return None

    def ancestors(self, go_id):
        """go_id and all its is_a ancestors in the full ontology."""
        if go_id not in self._ancestors:
            found = {go_id}
            for parent in self.terms[go_id]['is_a']:
                found |= self.ancestors(parent)
            self._ancestors[go_id] = frozenset(found)
        return self._ancestors[go_id]

    def label_matrix(self, pairs, row_of):
        """Propagated CSR uint8 labels from asserted (key, go_id) pairs: each term mapped (map_term), closed under is_a
        and restricted to the vocabulary. Returns (labels, the pairs whose term could not be mapped)."""
        rows, columns, dropped = [], [], []
        for key, go_id in pairs:
            term = self.map_term(go_id)
            if term is None:
                dropped.append((key, go_id))
                continue
            r = row_of[key]
            for ancestor in self.ancestors(term):
                j = self.index.get(ancestor)
                if j is not None:
                    rows.append(r)
                    columns.append(j)
        labels = sp.coo_matrix((np.ones(len(rows), np.int32), (rows, columns)), shape=(len(row_of), len(self.vocabulary))).tocsr()
        labels.data[:] = 1  # duplicates were summed
        return labels.astype(np.uint8), dropped

    def assert_closed(self, labels):
        """Fail if a positive label lacks a positive label on one of its in-vocabulary is_a parents."""
        labels = sp.csr_matrix(labels)
        for start in range(0, labels.shape[0], 4096):
            block = labels[start:start + 4096].toarray().astype(np.float32)
            gained = np.argwhere(self.propagate(block) != block)
            if len(gained):
                row, column = gained[0]
                raise ValueError(f'labels not closed under is_a: protein row {start + row} lacks '
                                 f'{self.vocabulary[column]}, an ancestor of one of its terms')

    def propagate(self, scores):
        """Every term gets the maximum score of itself and its descendants (is_a within the vocabulary)."""
        return self._propagation(np.asarray(scores, dtype=np.float32))


# ---------------------------------------------------------------------------------------------------------- fits
def fit_folder(root, method, training_set, seed):
    return Path(root) / 'runs' / method / training_set / f'seed{seed}'


def score_file(folder, split, side=None):
    """<split>_scores.npy, or <split>_scores_<side>.npy for a further checkpoint of the same fit (split: dev or test)."""
    return Path(folder) / f"{split}_scores{'_' + side if side else ''}.npy"


def write_fit(folder, dev_scores, test_scores, record, side=None):
    """Write a fit: raw scores (float32, rows in development/test order, columns in vocabulary order) of its primary
    checkpoint, and of further named checkpoints {name: (dev, test)}; then record.json, which marks the fit complete."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'record.json').unlink(missing_ok=True)
    for name, (dev, test) in {None: (dev_scores, test_scores), **(side or {})}.items():
        for split, x in (('dev', dev), ('test', test)):
            x = np.asarray(x, dtype=np.float32)
            if x.ndim != 2 or not np.isfinite(x).all():
                raise ValueError(f'{folder}: {split} scores {name or ""} of shape {x.shape} or not finite')
            np.save(score_file(folder, split, name), x)
    (folder / 'record.json').write_text(json.dumps(record, indent=1, default=str) + '\n')


def check_reproduced(saved, reloaded, tolerance=1e-5):
    """A reloaded checkpoint must give the development scores saved during training before it predicts the test set."""
    difference = float(np.abs(np.asarray(saved, np.float64) - np.asarray(reloaded, np.float64)).max())
    if not difference <= tolerance:  # NaN fails too
        raise ValueError(f'the reloaded checkpoint differs from its saved development scores by {difference}')
    return difference


# ---------------------------------------------------------------------------------------------------------- screen
def diamond():
    path = os.environ.get('DIAMOND') or shutil.which('diamond')
    if not path:
        raise SystemExit('DIAMOND not found: put `diamond` (version 2.2.4 in the paper) on PATH or set $DIAMOND')
    return path


def screen(root, proteins, work, threads=8):
    """Proteins ({id: sequence}) related to the benchmark: identical to a development or test sequence, or with a DIAMOND
    blastp hit to one in either direction (--very-sensitive, E <= 1e-3, no identity or coverage floor). Returns
    (exact ids, DIAMOND-hit ids)."""
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    held = {}
    for name in SPLITS.values():
        ds = load(root, name, labels=False)
        held.update({f'{name}:{h}': s for h, s in zip(ds.seq_ids, ds.sequences)})
    write_fasta(work / 'source.fasta', proteins.items())
    write_fasta(work / 'heldout.fasta', held.items())
    binary = diamond()
    for name in ('source', 'heldout'):
        subprocess.run([binary, 'makedb', '--in', work / f'{name}.fasta', '--db', work / name, '--threads', str(threads),
                        '--quiet'], check=True)
    flags = ['--very-sensitive', '--evalue', '1e-3', '--max-target-seqs', '0', '--dbsize', str(DIAMOND_DBSIZE),
             '--threads', str(threads), '--quiet', '--outfmt', '6', 'qseqid', 'sseqid', 'evalue', 'bitscore']
    for query, db in (('source', 'heldout'), ('heldout', 'source')):
        subprocess.run([binary, 'blastp', '--query', work / f'{query}.fasta', '--db', work / db,
                        '--out', work / f'{query}_vs_{db}.tsv', *flags], check=True)
    held_sequences = {s.upper() for s in held.values()}
    exact = {k for k, s in proteins.items() if s.upper() in held_sequences}
    with open(work / 'source_vs_heldout.tsv') as forward, open(work / 'heldout_vs_source.tsv') as reverse:
        hits = {line.split('\t')[0] for line in forward} | {line.split('\t')[1] for line in reverse}
    log(f'screen: {len(proteins):,} proteins, {len(exact):,} identical to a benchmark protein, {len(hits):,} with a '
        f'DIAMOND hit, {len(exact | hits):,} removed')
    return exact, hits
