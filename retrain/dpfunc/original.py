"""datasets/original_dpfunc: DPFunc's own molecular-function training data, rebuilt from the authors' release.

    python -m retrain.dpfunc.original [--swissprot FILE] [--threads 16]

Inputs:
- <data root>/external/dpfunc/dpfunc_published_dataset.zip, the authors' benchmark data (the Google Drive link in the
  README of github.com/CSUBioGroup/DPFunc): the archive of December 2025, checked by its SHA-256 (the authors have
  replaced this archive before; another one stops the build).
- --swissprot (default <data root>/external/uniprot/uniprot_sprot.fasta.gz): Swiss-Prot of UniProt release 2026_02,
  which the paper used to give the authors' entry names their accessions.

From the archive, LargeScaleData/mf_{train,test1,test2}.fasta,
*_go.txt and *_used_pid_list.pkl: the records the authors used (31,463 / 682 / 401, the paper's Table 3), as their
DataProcess/Process_data.ipynb selects them.

- The authors' validation (test1) and test (test2) proteins are training proteins here; nothing of theirs is held out.
- Proteins identical to, or with a DIAMOND hit to, a development or test protein are removed (retrain.common.screen).
- Records with the same sequence are one protein; its labels are the union of the records' labels.
- The supplied term sets are already propagated (on an older GO): every supplied term is taken as an assertion, mapped
  (alt_id to its primary id; an obsolete term with exactly one live replaced_by to it; other obsolete terms dropped),
  closed under is_a on the GO release of the data root and restricted to the vocabulary.
- accessions: the Swiss-Prot accession of each of the authors' entry names; for an entry renamed since, the Swiss-Prot
  accessions with exactly this sequence; else empty. The entry names are kept in the column entry_names.
"""
import argparse
import gzip
import hashlib
import io
import pickle
import zipfile
from collections import Counter
from pathlib import Path

from retrain.common import Ontology, add_argument, data_root, log, screen, seq_id, write_dataset
from retrain.dpfunc.common import external, inputs

SPLITS = (('train', 'train'), ('test1', 'valid'), ('test2', 'test'))  # archive name, the paper's name
ARCHIVE = 'dpfunc_published_dataset.zip'
ARCHIVE_SHA256 = '3fa11fc181b5db435c0ebe1ac6005ed38c1d0073283130b4b5767d163c1d56e2'


class _NoGlobals(pickle.Unpickler):
    """The used-pid lists are plain lists of strings; refuse anything that would import code."""

    def find_class(self, module, name):
        raise pickle.UnpicklingError(f'unexpected global {module}.{name} in a used-pid list')


def _fasta(text):
    out, name = {}, None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split()[0]
            if name in out:
                raise ValueError(f'duplicate FASTA id {name}')
            out[name] = []
        elif line:
            out[name].append(line)
    return {k: ''.join(v) for k, v in out.items()}


def records(root):
    """[(entry name, paper split, sequence, supplied GO terms)] of the records the authors used."""
    path = external(root) / ARCHIVE
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 24), b''):
            digest.update(block)
    if digest.hexdigest() != ARCHIVE_SHA256:
        raise SystemExit(f"{path} is not the authors' archive of December 2025 that the paper used (SHA-256 "
                         f'{digest.hexdigest()}, expected {ARCHIVE_SHA256})')
    out = []
    with zipfile.ZipFile(path) as z:
        for split, paper in SPLITS:
            sequences = _fasta(z.read(f'LargeScaleData/mf_{split}.fasta').decode())
            used = _NoGlobals(io.BytesIO(z.read(f'LargeScaleData/mf_{split}_used_pid_list.pkl'))).load()
            if not isinstance(used, list) or len(set(used)) != len(used) or not all(isinstance(p, str) for p in used):
                raise ValueError(f'{split}: the used-pid list is not a list of distinct strings')
            terms = {}
            for line in z.read(f'LargeScaleData/mf_{split}_go.txt').decode().splitlines():
                pid, go_id, aspect = line.split('\t')
                if aspect != 'mf':
                    raise ValueError(f'{split}: a non-MF line {line!r}')
                terms.setdefault(pid, set()).add(go_id)
            missing = [p for p in used if p not in sequences or p not in terms]
            if missing:
                raise ValueError(f'{split}: {len(missing)} used ids lack a sequence or labels, e.g. {missing[0]}')
            log(f'{split}: {len(used):,} records used of {len(sequences):,}')
            out += [(p, paper, sequences[p], sorted(terms[p])) for p in used]
    names = [r[0] for r in out]
    if len(set(names)) != len(names):
        raise ValueError('an entry name occurs in two splits')
    return out


def swissprot_accessions(path):
    """Swiss-Prot: entry name -> accession, and seq_id -> accessions (for entries renamed since the authors' release)."""
    by_name, by_sequence, accession, parts = {}, {}, None, []

    def flush():
        if accession:
            by_sequence.setdefault(seq_id(''.join(parts)), []).append(accession)
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                flush()
                _, accession, rest = line[1:].split('|', 2)
                by_name[rest.split()[0]] = accession
                parts = []
            else:
                parts.append(line.strip())
    flush()
    return by_name, by_sequence


def accessions_of(names, h, by_name, by_sequence):
    found = sorted({by_name[n] for n in names if n in by_name})
    return found or sorted(set(by_sequence.get(h, [])))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--swissprot', help='Swiss-Prot FASTA (gzip) of UniProt 2026_02 '
                                            '(default <data root>/external/uniprot/uniprot_sprot.fasta.gz)')
    parser.add_argument('--threads', type=int, default=16, help='DIAMOND threads')
    add_argument(parser)
    arguments = parser.parse_args()
    root = data_root(arguments.data)
    ontology = Ontology(root)
    recs = records(root)
    swissprot = arguments.swissprot or Path(root) / 'external' / 'uniprot' / 'uniprot_sprot.fasta.gz'
    by_name, by_sequence = swissprot_accessions(swissprot)
    unique = {}
    for _, _, sequence, _ in recs:
        unique.setdefault(seq_id(sequence), sequence)
    exact, hits = screen(root, unique, inputs(root) / 'original_screen', arguments.threads)
    removed = exact | hits

    proteins = {}  # seq_id -> record data, in the authors' order (train, then validation, then test)
    for name, split, sequence, terms in recs:
        h = seq_id(sequence)
        if h in removed:
            continue
        p = proteins.setdefault(h, dict(sequence=sequence, names=[], splits=[], terms=[]))
        p['names'].append(name)
        p['splits'].append(split)
        p['terms'].append((split, terms))
    ids = list(proteins)
    annotations, pairs, dropped = set(), [], Counter()
    for h in ids:
        for split, terms in proteins[h]['terms']:
            for go_id in terms:
                term = ontology.map_term(go_id)
                if term is None:
                    dropped[go_id] += 1
                    continue
                annotations.add((h, term, 'author_supplied', f'DPFunc mf_{split}'))
                pairs.append((h, term))
    labels, _ = ontology.label_matrix(pairs, {h: i for i, h in enumerate(ids)})
    log(f'{len(recs):,} records, {len(unique):,} distinct sequences, {len(removed):,} removed by the screen, '
        f'{len(ids):,} proteins; {sum(dropped.values()):,} supplied terms dropped (obsolete without a single replacement)')
    accessions = [accessions_of(proteins[h]['names'], h, by_name, by_sequence) for h in ids]
    log(f'accessions: {sum(map(bool, accessions)):,} of {len(ids):,} proteins')
    write_dataset(root, 'original_dpfunc', ids, accessions, [proteins[h]['sequence'] for h in ids], annotations,
                  labels, ontology, columns=dict(entry_names=[','.join(proteins[h]['names']) for h in ids],
                                                 source_splits=[','.join(proteins[h]['splits']) for h in ids]))


if __name__ == '__main__':
    main()
