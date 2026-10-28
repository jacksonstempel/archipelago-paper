"""DPFunc's model inputs, made by the same code for every training set, the development set and the test set.

    python -m retrain.dpfunc.inputs graphs [--datasets ...] [--cpus 16]
    python -m retrain.dpfunc.inputs bank {original,archipelago,evaluation} [--cpus 16]
    python -m retrain.dpfunc.inputs step DATASET

graphs  contact graph of each protein from its backbone in the structure bank (retrain.structures): directed edges
        between C-alpha atoms at most 10 Angstrom apart (the paper; the released code's default is 12), no self loops,
        ordered by (source, target) as the released double loop makes them -> inputs/dpfunc/graphs/<ab>/<seq_id>.npy.
bank    one memory-mapped file of residue features and one of edges for a group of datasets (original: original_dpfunc;
        archipelago, which also holds goa, hq and the other HQ subsets; evaluation: development then test), in dataset
        order: banks/<bank>/{seq_ids.json, feature_offsets.npy, features.npy, edge_offsets.npy, edges.npy}. Above 2,700
        residues the graph and the features cover residues 1-1,400 (the structure bank's rule). ArchipelaGO's bank is
        large (about 0.6 TB of FP32 features).
step    a dataset's rows in its bank and its InterPro count matrix on the released 22,369-entry vocabulary (entries
        outside it are dropped, as the released model does) -> steps/<dataset>/{rows.npy, domains.npz}.

Needs inputs/dpfunc/esm1b (retrain.dpfunc.embed), inputs/dpfunc/domains/counts.json.gz (retrain.dpfunc.domains) and
external/dpfunc/inter_idx.pkl (data/inter_idx.pkl of github.com/CSUBioGroup/DPFunc at commit c5faf4df).
"""
import argparse
import io
import json
import multiprocessing
import pickle
import shutil
from functools import partial

import numpy as np
import scipy.sparse as sp

from retrain.common import add_argument, data_root, load, log
from retrain.dpfunc.common import (BANKS, CONTACT_ANGSTROM, DIMENSION, N_INTERPRO, bank_of, esm_file, external, inputs,
                                   node_positions)
from retrain.dpfunc.domains import DATASETS, read_counts
from retrain.structures.assemble import backbone_path
from retrain.structures.common import PROXIES


def edges(ca, threshold=CONTACT_ANGSTROM):
    """Directed C-alpha contacts at most `threshold` apart, no self loops, ordered by (source, target): the released
    double loop over residue pairs."""
    from scipy.spatial import cKDTree
    ca = np.asarray(ca, dtype=np.float64)
    pairs = cKDTree(ca).query_pairs(threshold + 1e-9, output_type='ndarray')
    if len(pairs):
        delta = ca[pairs[:, 0]] - ca[pairs[:, 1]]
        pairs = pairs[np.sqrt(np.sum(delta * delta, axis=1)) <= threshold]
    pairs = np.concatenate([pairs, pairs[:, ::-1]], axis=0)
    if len(pairs):
        pairs = pairs[np.lexsort((pairs[:, 1], pairs[:, 0]))]
    return pairs.T.astype(np.int64).reshape(2, -1)


def backbone_ca(root, seq_id, sequence):
    """C-alpha coordinates of the protein's structure; the structure must cover exactly node_positions(length) and
    match the protein's residues (up to the stated proxies)."""
    with np.load(backbone_path(root, seq_id), allow_pickle=False) as z:
        ca, stored, full_length = np.asarray(z['ca'], dtype=np.float32), str(z['sequence']), int(z['full_length'])
    if full_length != len(sequence) or len(ca) != len(node_positions(len(sequence))) or len(stored) != len(ca):
        raise ValueError(f'{seq_id}: the structure does not cover the residues DPFunc uses')
    if any(s != r and s != PROXIES.get(r) for s, r in zip(stored, sequence)):
        raise ValueError(f'{seq_id}: structure residues differ from the protein beyond the stated proxies')
    return ca


def graph_file(root, seq_id):
    return inputs(root) / 'graphs' / seq_id[:2] / f'{seq_id}.npy'


def make_graph(root, item):
    seq_id, sequence = item
    target = graph_file(root, seq_id)
    if not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = target.with_name(f'{seq_id}.tmp.npy')
        np.save(temporary, edges(backbone_ca(root, seq_id, sequence)))
        temporary.replace(target)


def graphs(root, arguments):
    sequences = {}
    for name in arguments.datasets:
        ds = load(root, name, labels=False)
        sequences.update(zip(ds.seq_ids, ds.sequences))
    todo = [(h, s) for h, s in sorted(sequences.items()) if not graph_file(root, h).exists()]
    log(f'graphs: {len(todo):,} of {len(sequences):,} proteins to make')
    with multiprocessing.Pool(arguments.cpus) as pool:
        for i, _ in enumerate(pool.imap_unordered(partial(make_graph, root), todo, chunksize=64), 1):
            if i % 20000 == 0 or i == len(todo):
                log(f'graphs: {i:,}/{len(todo):,}')


def _fill_features(job):
    root, path, seq_ids, lengths, offsets = job
    features = np.load(path, mmap_mode='r+')
    for h, n, a, b in zip(seq_ids, lengths, offsets[:-1], offsets[1:]):
        x = np.load(esm_file(root, h), mmap_mode='r')
        if x.dtype != np.float32 or x.shape != (n, DIMENSION):
            raise ValueError(f'ESM-1b features of {h}: {x.dtype} {x.shape}, sequence length {n}')
        features[a:b] = x[node_positions(n)]
    features.flush()
    return len(seq_ids)


def _fill_edges(job):
    root, path, seq_ids, offsets = job
    out = np.load(path, mmap_mode='r+')
    for h, a, b in zip(seq_ids, offsets[:-1], offsets[1:]):
        out[:, a:b] = np.load(graph_file(root, h))
    out.flush()
    return len(seq_ids)


def bank(root, arguments):
    name = arguments.name
    seq_ids, lengths = [], []
    for dataset in BANKS[name]:
        ds = load(root, dataset, labels=False)
        seq_ids += ds.seq_ids
        lengths += [len(s) for s in ds.sequences]
    if len(set(seq_ids)) != len(seq_ids):
        raise ValueError(f'bank {name}: a protein repeats across {BANKS[name]}')
    missing = [h for h in seq_ids if not esm_file(root, h).exists() or not graph_file(root, h).exists()]
    if missing:
        raise SystemExit(f'bank {name}: {len(missing):,} proteins lack ESM-1b features or a graph, e.g. {missing[0]}')
    folder = inputs(root) / 'banks' / name
    building = folder.with_name(f'.{name}.building')
    shutil.rmtree(building, ignore_errors=True)
    building.mkdir(parents=True)
    (building / 'seq_ids.json').write_text(json.dumps(seq_ids))
    offsets = np.r_[0, np.cumsum([len(node_positions(n)) for n in lengths])].astype(np.int64)
    np.save(building / 'feature_offsets.npy', offsets)
    log(f'bank {name}: {len(seq_ids):,} proteins, {int(offsets[-1]):,} residues, '
        f'{int(offsets[-1]) * DIMENSION * 4 / 2 ** 30:,.0f} GiB of features')
    path = building / 'features.npy'
    np.lib.format.open_memmap(path, mode='w+', dtype=np.float32, shape=(int(offsets[-1]), DIMENSION)).flush()
    jobs = [(root, path, seq_ids[a:a + 1000], lengths[a:a + 1000], offsets[a:a + 1001]) for a in range(0, len(seq_ids), 1000)]
    done = 0
    with multiprocessing.Pool(arguments.cpus) as pool:
        for n in pool.imap_unordered(_fill_features, jobs):
            done += n
            log(f'bank {name}: features of {done:,}/{len(seq_ids):,} proteins')
    counts = [np.load(graph_file(root, h), mmap_mode='r').shape[1] for h in seq_ids]
    edge_offsets = np.r_[0, np.cumsum(counts)].astype(np.int64)
    np.save(building / 'edge_offsets.npy', edge_offsets)
    path = building / 'edges.npy'
    np.lib.format.open_memmap(path, mode='w+', dtype=np.int64, shape=(2, int(edge_offsets[-1]))).flush()
    jobs = [(root, path, seq_ids[a:a + 2000], edge_offsets[a:a + 2001]) for a in range(0, len(seq_ids), 2000)]
    with multiprocessing.Pool(arguments.cpus) as pool:
        sum(pool.imap_unordered(_fill_edges, jobs))
    shutil.rmtree(folder, ignore_errors=True)
    building.rename(folder)
    log(f'bank {name}: {int(edge_offsets[-1]):,} edges written')


class _NoGlobals(pickle.Unpickler):
    def find_class(self, module, name):
        raise pickle.UnpicklingError(f'unexpected global {module}.{name}')


def interpro_vocabulary(root):
    """{InterPro entry: column} of the released model (data/inter_idx.pkl of the DPFunc repository)."""
    index = _NoGlobals(io.BytesIO((external(root) / 'inter_idx.pkl').read_bytes())).load()
    if not isinstance(index, dict):
        index = {entry: i for i, entry in enumerate(index)}
    if len(index) != N_INTERPRO or set(index.values()) != set(range(N_INTERPRO)):
        raise ValueError('inter_idx.pkl is not the released 22,369-entry InterPro vocabulary')
    return index


def step(root, arguments):
    name = arguments.name
    ds = load(root, name, labels=False)
    bank_ids = json.loads((inputs(root) / 'banks' / bank_of(name) / 'seq_ids.json').read_text())
    row_of = {h: i for i, h in enumerate(bank_ids)}
    missing = [h for h in ds.seq_ids if h not in row_of]
    if missing:
        raise ValueError(f'{name}: {len(missing):,} proteins are not in bank {bank_of(name)}')
    counts, vocabulary = read_counts(root), interpro_vocabulary(root)
    indices, values, indptr, pairs, dropped = [], [], [0], 0, 0
    for h in ds.seq_ids:
        kept = sorted((vocabulary[t], n) for t, n in counts[h].items() if t in vocabulary)
        pairs += len(counts[h])
        dropped += len(counts[h]) - len(kept)
        indices += [j for j, _ in kept]
        values += [n for _, n in kept]
        indptr.append(len(indices))
    matrix = sp.csr_matrix((np.asarray(values, np.float32), np.asarray(indices, np.int64), np.asarray(indptr, np.int64)),
                           shape=(len(ds.seq_ids), N_INTERPRO))
    folder = inputs(root) / 'steps' / name
    folder.mkdir(parents=True, exist_ok=True)
    np.save(folder / 'rows.npy', np.asarray([row_of[h] for h in ds.seq_ids], dtype=np.int64))
    sp.save_npz(folder / 'domains.npz', matrix, compressed=True)
    log(f'step {name}: {len(ds.seq_ids):,} proteins in bank {bank_of(name)}; {dropped:,} of {pairs:,} protein-entry pairs '
        f'({dropped / max(pairs, 1):.1%}) name entries outside the released InterPro vocabulary and are dropped')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('stage', choices=['graphs', 'bank', 'step'])
    parser.add_argument('name', nargs='?', help='bank: the bank; step: the dataset')
    parser.add_argument('--datasets', nargs='+', default=list(DATASETS), help='graphs: the datasets to cover')
    parser.add_argument('--cpus', type=int, default=16)
    add_argument(parser)
    arguments = parser.parse_args()
    if arguments.stage == 'bank' and arguments.name not in BANKS:
        parser.error(f'bank needs one of {list(BANKS)}')
    if arguments.stage == 'step' and not arguments.name:
        parser.error('step needs a dataset')
    globals()[arguments.stage](data_root(arguments.data), arguments)


if __name__ == '__main__':
    main()
