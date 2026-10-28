"""Shared by the structure-bank stages: the bank's folder, the modeled sequence of a protein, PDB backbones, small tables.

The bank lives in <data root>/inputs/structures/. Every protein is keyed by its seq_id (SHA-256 of the upper-case sequence).
"""
import gzip
import json
import os
import time
from pathlib import Path

import numpy as np

LONG_PROTEIN = 2700  # above this length only residues 1-1,400 are modeled (AlphaFold DB's F1 fragment of long proteins)
FRAGMENT = 1400
STANDARD = set('ACDEFGHIKLMNPQRSTVWY')
PROXIES = dict(U='C', O='K', B='D', Z='E', X='A')  # stand-ins for residues AlphaFold2 has no atom type for; positions kept
THREE = dict(ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLN='Q', GLU='E', GLY='G', HIS='H', ILE='I', LEU='L', LYS='K',
             MET='M', PHE='F', PRO='P', SER='S', THR='T', TRP='W', TYR='Y', VAL='V', SEC='U', PYL='O', ASX='B', GLX='Z')


def bank(root):
    return Path(root) / 'inputs' / 'structures'


def modeled_raw(sequence):
    """(the residues a structure must cover, coverage): the whole sequence, or residues 1-1,400 above 2,700 residues."""
    if len(sequence) > LONG_PROTEIN:
        return sequence[:FRAGMENT], 'residues_1_1400'
    return sequence, 'full_chain'


def target_sequence(sequence):
    """The sequence ColabFold folds: modeled_raw() with U/O/B/Z/X replaced by their proxies."""
    raw, _ = modeled_raw(sequence)
    return ''.join(a if a in STANDARD else PROXIES.get(a, 'A') for a in raw)


def read_structure(payload):
    """(residue sequence, N/CA/C coordinates (L, 3, 3), per-residue B-factor = pLDDT) of a one-chain PDB file (bytes, gzip
    or plain). Fails on a second chain, a residue without N, CA or C, insertion codes or non-finite values."""
    if payload[:2] == b'\x1f\x8b':
        payload = gzip.decompress(payload)
    residues, chains, order = {}, set(), []
    for line in payload.decode().splitlines():
        if line.startswith('ENDMDL'):
            break
        if not line.startswith('ATOM'):
            continue
        chains.add(line[21])
        key = (line[22:27], line[17:20])
        if key not in residues:
            residues[key] = {}
            order.append(key)
        atom = line[12:16].strip()
        if atom in ('N', 'CA', 'C'):
            residues[key][atom] = (float(line[30:38]), float(line[38:46]), float(line[46:54]), float(line[60:66]))
    if len(chains) != 1:
        raise ValueError(f'expected one chain, found {sorted(chains)}')
    if any(k[0][4] != ' ' for k in order):
        raise ValueError('insertion codes present')
    missing = [i for i, k in enumerate(order) if len(residues[k]) != 3]
    if missing:
        raise ValueError(f'{len(missing)} residues lack N, CA or C (first at index {missing[0]})')
    sequence = ''.join(THREE.get(k[1], 'X') for k in order)
    xyz = np.array([[residues[k][a][:3] for a in ('N', 'CA', 'C')] for k in order], dtype=np.float32)
    plddt = np.array([residues[k]['CA'][3] for k in order], dtype=np.float32)
    if not (np.isfinite(xyz).all() and np.isfinite(plddt).all()):
        raise ValueError('non-finite coordinates')
    return sequence, xyz, plddt


def read_tsv(path):
    with open(path) as handle:
        header = handle.readline().rstrip('\n').split('\t')
        return [dict(zip(header, line.rstrip('\n').split('\t'))) for line in handle if line.strip()]


def write_tsv(path, rows, columns):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f'{path.name}.{os.getpid()}.tmp')
    with open(temporary, 'w') as handle:
        handle.write('\t'.join(columns) + '\n')
        handle.writelines('\t'.join(str(r.get(c, '')) for c in columns) + '\n' for r in rows)
    temporary.replace(path)


def write_bytes(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f'{path.name}.{os.getpid()}.tmp')
    temporary.write_bytes(payload)
    temporary.replace(path)


def read_jsonl(path):
    return [json.loads(line) for line in open(path)] if Path(path).exists() else []


class Progress:
    """Rate and ETA lines for long loops (at most one per `every` seconds)."""

    def __init__(self, label, total, every=60):
        self.label, self.total, self.every = label, total, every
        self.started = self.last = time.monotonic()
        self.done = 0

    def step(self, n=1, **extra):
        from retrain.common import log
        self.done += n
        now = time.monotonic()
        if now - self.last >= self.every or self.done == self.total:
            self.last = now
            rate = self.done / max(now - self.started, 1e-9)
            eta = (self.total - self.done) / rate / 60 if rate else float('nan')
            log(f'{self.label}: {self.done:,}/{self.total:,}, {rate:.2f}/s, ETA {eta:.1f} min', json.dumps(extra) if extra else '')
