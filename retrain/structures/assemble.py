"""Structure bank, part 3: one backbone file per protein and the bank's index.

    python -m retrain.structures.assemble [--workers 8]

For every protein of requests.tsv: its AlphaFold DB model (afdb/selected.tsv) or else its ColabFold fold, checked residue by
residue against the protein (AlphaFold DB: identical to the modeled residues; ColabFold: identical up to the stated
U/O/B/Z/X proxies), written as backbone/<seq_id>.npz (n, ca, c: float32 (L, 3); sequence: the structure's residues;
plddt; full_length: the protein's length). Residue i of a backbone is residue i + 1 of the protein; above 2,700 residues
the backbone covers residues 1-1,400. index.tsv lists every protein's source; it is written only when every protein has
a structure.
"""
import argparse
import hashlib
from collections import Counter
from multiprocessing import Pool

import numpy as np

from retrain.common import add_argument, data_root, log
from retrain.structures.afdb import load_requests
from retrain.structures.colabfold import model_path
from retrain.structures.common import (PROXIES, Progress, bank, modeled_raw, read_structure, read_tsv, target_sequence,
                                       write_tsv)

COLUMNS = ['seq_id', 'length', 'sets', 'source', 'model', 'coverage', 'mean_plddt']


def backbone_path(root, seq_id):
    return bank(root) / 'backbone' / f'{seq_id}.npz'


def one(job):
    root, r, model_row = job
    sequence = r['sequence']
    raw, coverage = modeled_raw(sequence)
    if model_row is not None:
        source, model = 'afdb', f"{model_row['entry']}-model_v{model_row['version']}"
        path, expected = bank(root) / model_row['path'], raw
    else:
        target = target_sequence(sequence)
        key = hashlib.sha256(target.encode()).hexdigest()
        source, model, path, expected = 'colabfold', key, model_path(root, key), target
        if not path.exists():
            return None
    observed, xyz, plddt = read_structure(path.read_bytes())
    if observed != expected:
        raise ValueError(f"{r['seq_id']}: structure residues differ from the modeled sequence ({source} {model})")
    if any(o != a and o != PROXIES.get(a) for o, a in zip(observed, raw)):
        raise ValueError(f"{r['seq_id']}: structure residues differ from the protein beyond the proxies")
    out = backbone_path(root, r['seq_id'])
    temporary = out.with_name(f'{out.name}.tmp.npz')
    np.savez_compressed(temporary, n=xyz[:, 0], ca=xyz[:, 1], c=xyz[:, 2], sequence=np.array(observed), plddt=plddt,
                        full_length=np.int64(len(sequence)))
    temporary.replace(out)
    return dict(seq_id=r['seq_id'], length=len(sequence), sets=r['sets'], source=source, model=model, coverage=coverage,
                mean_plddt=round(float(plddt.mean()), 2))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--workers', type=int, default=8)
    add_argument(parser)
    arguments = parser.parse_args()
    root = data_root(arguments.data)
    (bank(root) / 'backbone').mkdir(parents=True, exist_ok=True)
    selected = {r['seq_id']: r for r in read_tsv(bank(root) / 'afdb' / 'selected.tsv')}
    jobs = [(root, r, selected.get(r['seq_id'])) for r in load_requests(root)]
    rows, missing, progress = [], [], Progress('assemble', len(jobs))
    with Pool(arguments.workers) as pool:
        for job, row in zip(jobs, pool.imap(one, jobs, chunksize=64)):
            if row is None:
                missing.append(job[1]['seq_id'])
            else:
                rows.append(row)
            progress.step()
    if missing:
        raise SystemExit(f'{len(missing):,} proteins have no structure yet (first {missing[0]}); run the ColabFold stages')
    write_tsv(bank(root) / 'index.tsv', rows, COLUMNS)
    log('assemble:', dict(Counter(r['source'] for r in rows)), 'residues 1-1,400 only:',
        sum(r['coverage'] != 'full_chain' for r in rows))


if __name__ == '__main__':
    main()
