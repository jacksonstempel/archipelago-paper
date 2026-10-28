"""datasets/original_transfew: TransFew's own MF training data, rebuilt from the authors' release.

Input: the authors' data release TFewData.tar.gz (https://doi.org/10.5281/zenodo.22139108; 20,990,155,097 bytes,
SHA-256 2f9007f96b1fdd28cb12c35e452f4415eaa69e8e6b2ae4e5026d100de672bf40), extracted so that
<root>/external/transfew/TFewData/ exists. Files read:
  mf/train_data.pickle, mf/validation_data.pickle   the training and validation records: keys 'protein' (accessions)
                                                     and 'labels' (one binary (1, L) tensor per protein over
                                                     mf/sorted_terms, the authors' terms, already propagated); record
                                                     ids 'train:<position>:<accession>', in file order
  mf/{sorted_terms,term_indicies,train_proteins,validation_proteins,all_proteins}.pickle   index files (cross-check)
  train_sequences.fasta, test_fasta.fasta           the authors' sequences of their training and test proteins
  test/output_t1_t2/groundtruths/full/mf_{swissprot,trembl}.tsv   their MF test proteins' terms

The dataset:
  - Records: the authors' train (70,773) and validation (7,864) records with their released terms, plus their MF test
    proteins (Swiss-Prot and TrEMBL) with their ground-truth terms. The authors' validation and test proteins are
    training data here; evaluation uses the paper's benchmark.
  - Sequences: the authors' own FASTA files, looked up by the record's accession (exactly one sequence per accession,
    or the build fails). No Swiss-Prot release or accession aliases are needed: the paper's build also took every
    sequence from these files.
  - The test_data.pickle of the release has no labels and is not used; the 78,637 train and validation records are
    checked by count and the three defining pickles by SHA-256.
  - Screen: a record is removed if its sequence equals a development or test sequence or has a DIAMOND hit to one
    (retrain.common.screen).
  - One row per distinct sequence: records with the same sequence are merged, their terms united.
  - Labels: every released term mapped on the paper's GO release (alt_id -> primary id; an obsolete term with exactly
    one live replaced_by -> that term; other obsolete terms dropped), closed under is_a, restricted to the vocabulary.
    Proteins left without a vocabulary label are dropped.
Expected result (the paper's set): 70,119 proteins and 533,182 propagated label pairs. The paper's build read the
same records through an intermediate extraction of these pickles in their file order; if the counts match but a
retrained model differs slightly, row order (and with it the training shuffle) is the first thing to compare.

    python -m retrain.transfew.original [--data DIR] [--threads 16]

CPU only. The two record pickles also hold the authors' embeddings (16.8 GB together) and are unpickled whole, so the
build needs about 40 GB of memory. They contain PyTorch storages, so they are read with pickle: open only the release
files from the DOI above.
"""
import argparse
import collections
import hashlib
import io
import pickle
from pathlib import Path

import numpy as np

from retrain import common

SPLITS = ('train', 'validation')
TEST_PARTS = ('swissprot', 'trembl')
RELEASED_RECORDS = 78_637  # 70,773 train + 7,864 validation
CHECKSUMS = {  # SHA-256 of the release members that define the records
    'mf/sorted_terms.pickle': 'ebb8b3a8517b537f95b8b3c9b4947b2b768284fa0c6e420d047cfe75af9b8afe',
    'mf/train_data.pickle': 'c9af1f8874e386d79ff6e9a3997f0168c8906a0af6579f1c455d5ebf6b4f4c86',
    'mf/validation_data.pickle': '10ca41bd176b1f8c3c409e255770a9e28c6c4ab42e37edcdd7be0d91bc9c188b',
}


class _IndexUnpickler(pickle.Unpickler):
    """The index pickles hold only builtin containers, strings and numbers; anything else is refused."""

    def find_class(self, module, name):
        if module == 'builtins' and name in ('set', 'frozenset', 'list', 'dict', 'tuple'):
            return getattr(__import__('builtins'), name)
        raise pickle.UnpicklingError(f'refusing {module}.{name}')


class _RecordUnpickler(pickle.Unpickler):
    """The record pickles hold PyTorch tensors saved on the authors' GPUs; load them on the CPU."""

    def find_class(self, module, name):
        if module == 'torch.storage' and name == '_load_from_bytes':
            import torch
            return lambda b: torch.load(io.BytesIO(b), map_location='cpu', weights_only=False)
        return super().find_class(module, name)


def index_pickle(release, name):
    with open(release / 'mf' / f'{name}.pickle', 'rb') as handle:
        return _IndexUnpickler(handle).load()


def check_files(release):
    """The release files whose content defines the records must be the ones the paper used."""
    for name, expected in CHECKSUMS.items():
        digest = hashlib.sha256()
        with open(release / name, 'rb') as handle:
            for block in iter(lambda: handle.read(1 << 24), b''):
                digest.update(block)
        if digest.hexdigest() != expected:
            raise ValueError(f'{release / name}: SHA-256 {digest.hexdigest()}, expected {expected}')


def released_records(release):
    """[(split, protein id, released terms)] of the train and validation pickles, in file order."""
    check_files(release)
    terms = index_pickle(release, 'sorted_terms')
    out = []
    for split in SPLITS:
        common.log(f'reading mf/{split}_data.pickle')
        with open(release / 'mf' / f'{split}_data.pickle', 'rb') as handle:
            data = _RecordUnpickler(handle).load()
        proteins, labels = data['protein'], data['labels']
        if len(proteins) != len(labels):
            raise ValueError(f'{split}_data: {len(proteins)} proteins but {len(labels)} label vectors')
        for protein, vector in zip(proteins, labels):
            vector = np.asarray(vector).reshape(-1)
            if len(vector) != len(terms) or not np.isin(vector, (0, 1)).all():
                raise ValueError(f'{split}_data {protein}: label vector of length {len(vector)} or not binary')
            out.append((split, protein, [terms[j] for j in np.flatnonzero(vector)]))
        del data
    if len(out) != RELEASED_RECORDS:
        raise ValueError(f'{len(out):,} train and validation records, expected {RELEASED_RECORDS:,}')
    check_index(release, out, terms)
    return out


def check_index(release, records, terms):
    """The records must match the release's own index files: split membership, and every term's protein count against
    term_indicies (Preprocess.py: count > t and < 70,000 for t in 0, 30, 50, 100, 250, positions in sorted_terms).
    A mismatch means the pickles were misread."""
    by_split, counts = collections.defaultdict(set), collections.Counter()
    for split, protein, released in records:
        by_split[split].add(protein)
        counts.update(released)
    position = {t: i for i, t in enumerate(terms)}
    indices = index_pickle(release, 'term_indicies')
    problems = [f'{s} proteins' for s in SPLITS if by_split[s] != set(index_pickle(release, f'{s}_proteins'))]
    if by_split['train'] | by_split['validation'] != set(index_pickle(release, 'all_proteins')):
        problems.append('all proteins')
    for t in (0, 30, 50, 100, 250):
        if sorted(position[g] for g, c in counts.items() if t < c < 70000) != sorted(indices[t]):
            problems.append(f'term_indicies[{t}]')
    if problems:
        raise ValueError(f'released records disagree with the release index files: {problems}')


def fasta_by_accession(path):
    """{accession: set of sequences}; headers 'sp|P12345|NAME ...' or 'P12345 ...'. An accession may repeat (its
    sequences are collected; one_sequence() decides)."""
    out, header, parts = collections.defaultdict(set), None, []

    def flush():
        if header is not None:
            first = header.split()[0]
            out[first.split('|')[1] if first.count('|') >= 2 else first].add(''.join(parts))

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                flush()
                header, parts = line[1:], []
            elif line:
                parts.append(line)
    flush()
    return out


def one_sequence(sequences, accession, what):
    found = sequences.get(accession, set())
    if len(found) != 1:
        raise ValueError(f'{what} {accession}: {len(found)} sequences in the authors\' FASTA')
    return next(iter(found))


def records(release):
    """Every record (training, validation, test) as dict(record, accession, split, sequence, terms)."""
    out = []
    authors = fasta_by_accession(release / 'train_sequences.fasta')
    position = collections.Counter()
    for split, protein, released in released_records(release):
        out.append(dict(record=f'{split}:{position[split]}:{protein}', accession=protein, split=split,
                        sequence=one_sequence(authors, protein, f'{split} record'), terms=released))
        position[split] += 1
    test_sequences = fasta_by_accession(release / 'test_fasta.fasta')
    for part in TEST_PARTS:
        terms = collections.defaultdict(list)
        with open(release / f'test/output_t1_t2/groundtruths/full/mf_{part}.tsv') as handle:
            for line in handle:
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 2 or not fields[1].startswith('GO:'):
                    raise ValueError(f'mf_{part}.tsv: unexpected line {line!r}')
                terms[fields[0]].append(fields[1])
        for i, (protein, ts) in enumerate(terms.items()):
            out.append(dict(record=f'test_{part}:{i}:{protein}', accession=protein, split=f'test_{part}',
                            sequence=one_sequence(test_sequences, protein, f'{part} test protein'),
                            terms=sorted(set(ts))))
    common.log('records:', dict(collections.Counter(r['split'] for r in out)))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    common.add_argument(parser)
    parser.add_argument('--threads', type=int, default=16, help='DIAMOND threads')
    args = parser.parse_args()
    root = common.data_root(args.data)
    release = root / 'external' / 'transfew' / 'TFewData'
    if (root / 'datasets' / 'original_transfew').exists():
        raise FileExistsError(root / 'datasets' / 'original_transfew')
    rows = records(release)
    unique = {}
    for r in rows:
        unique.setdefault(common.seq_id(r['sequence']), r['sequence'])
    exact, hits = common.screen(root, unique, root / 'inputs' / 'transfew' / 'original_screen', args.threads)
    removed = exact | hits
    kept = [r for r in rows if common.seq_id(r['sequence']) not in removed]

    # one row per distinct sequence, in record order (train, validation, test); terms united over its records
    members = collections.defaultdict(list)
    for r in kept:
        members[common.seq_id(r['sequence'])].append(r)
    order = list(members)
    ontology = common.Ontology(root)
    pairs = [(h, t) for h in order for t in sorted({t for r in members[h] for t in r['terms']})]
    labels, dropped = ontology.label_matrix(pairs, {h: i for i, h in enumerate(order)})
    keep = np.flatnonzero(labels.getnnz(axis=1) > 0)
    order, labels = [order[i] for i in keep], labels[keep]
    annotations = []
    for h in order:
        for r in members[h]:
            for t in r['terms']:
                term = ontology.map_term(t)
                if term is not None:
                    source = f"TransFew release {r['split']} {r['accession']}" + (f' (released as {t})' if term != t else '')
                    annotations.append((h, term, 'released, already propagated by the authors', source))
    common.write_dataset(
        root, 'original_transfew', order, [sorted({r['accession'] for r in members[h]}) for h in order],
        [members[h][0]['sequence'] for h in order], annotations, labels, ontology,
        columns=dict(records=[','.join(r['record'] for r in members[h]) for h in order]))
    common.log(f"original_transfew: {len(rows):,} records, {len(rows) - len(kept):,} removed by the screen, "
               f"{len(order):,} proteins ({len(members) - len(order)} without a vocabulary label dropped), "
               f"{len({t for _, t in dropped}):,} obsolete terms without a single replacement dropped "
               f"({len(dropped):,} protein-term pairs)")


if __name__ == '__main__':
    main()
