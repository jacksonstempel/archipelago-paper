"""hq_random: HQ plus 150,000 proteins drawn uniformly at random from all computationally annotated TrEMBL proteins.

    python -m retrain.datasets.random_draw [--data DIR]                      the released draw (default)
    python -m retrain.datasets.random_draw --from-pool [--threads N] [--data DIR]   draw again from the pool

Default: the 150,000 drawn proteins as released, under <data root>/release/hq_random/, in the format of the candidates'
release files (retrain.datasets.candidates): proteins.parquet (protein_id, accessions, sequence, go_terms; rows in draw
order) and annotations.parquet (protein_id, go_id, evidence, assigned_by). The same checks apply (each sequence hashes
to its protein_id, the annotations propagate to exactly go_terms). The rank column order_rank is the position in that
list.

How the list was drawn (--from-pool repeats it). The pool: every TrEMBL sequence with molecular-function annotations and
no experimental or TAS/IC support (77,896,027 distinct sequences in the paper's GOA release), not the cluster
representatives: random representatives would inherit the diversity that clustering selects for. HQ, development and
test sequences are removed first. The pool's seq_ids are sorted, then permuted with numpy default_rng(20261023); the
first 300,000 drawn are checked in draw order, and the first 150,000 that pass two filters are kept:
  screen      no exact or DIAMOND relation to a development or test protein (retrain.common.screen)
  vocabulary  every is_a ancestor of every annotated term lies in the output vocabulary (as for the candidates)
No other selection is made. The rank column order_rank is then the position in the draw, skipped draws included.
Its input, under <data root>/external/computational_pool/ (not released; about 20 GB):
  proteins.parquet    seq_id, sequence, accessions (comma-separated): one row per distinct sequence
  assertions.parquet  seq_id, go_id, evidence, assigned_by: the computational MF annotations of each sequence
Reading the pool's 77.9 million identifiers needs about 16 GB of memory; the screen of 300,000 sequences takes minutes
with 16 threads.
"""
import argparse
from collections import defaultdict

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as pads
import pyarrow.parquet as pq

from retrain.common import Ontology, add_argument, data_root, load, log, screen, seq_id
from retrain.datasets.added import write_added
from retrain.datasets.candidates import annotation_row, checked

SEED, DRAW, KEEP = 20261023, 300_000, 150_000
POOL = 77_896_027  # distinct sequences of the paper's computational pool
ROOT_TERM = 'GO:0003674'


def from_release(root, ontology):
    """The released draw: (protein ids, ranks, accessions, sequences, annotation rows, labels)."""
    source = root / 'release' / 'hq_random'
    labels, ids, annotations = checked(source, ontology)
    table = pq.read_table(source / 'proteins.parquet', columns=['protein_id', 'accessions', 'sequence']).to_pydict()
    if len(table['protein_id']) != KEEP:
        raise ValueError(f"{source}: {len(table['protein_id']):,} proteins; the paper drew {KEEP:,}")
    wrong = [h for h, s in zip(table['protein_id'], table['sequence']) if seq_id(s) != h]
    if wrong:
        raise ValueError(f'{len(wrong)} sequences do not hash to their protein_id, e.g. {wrong[0]}')
    columns = [annotations.column(k).to_pylist() for k in ('protein_id', 'go_id', 'evidence', 'assigned_by')]
    return (table['protein_id'], range(1, KEEP + 1), [a or [] for a in table['accessions']],
            [s.upper() for s in table['sequence']], [annotation_row(*r) for r in zip(*columns)], labels)


def fixed_width_ids(chunked):
    """The 64-character hex seq_ids of a string column as one numpy 'S64' array (no Python strings). Chunk by chunk: the
    whole column (5 GB) would overflow a string array's 32-bit offsets."""
    parts = []
    for array in chunked.chunks:
        large = pa.types.is_large_string(array.type)
        offsets = np.frombuffer(array.buffers()[1], dtype=np.int64 if large else np.int32, count=len(array) + 1,
                                offset=array.offset * (8 if large else 4))
        if array.null_count or not (np.diff(offsets) == 64).all():
            raise ValueError('pool seq_ids are not 64-character SHA-256 hex digests')
        parts.append(np.frombuffer(array.buffers()[2], dtype='S64', count=len(array), offset=int(offsets[0])))
    return np.concatenate(parts)


def draw_order(pool, excluded):
    """The first DRAW seq_ids of the draw: sorted pool minus the excluded ids, permuted."""
    ids = fixed_width_ids(pq.read_table(pool / 'proteins.parquet', columns=['seq_id']).column('seq_id'))
    if len(np.unique(ids)) != len(ids) or len(ids) != POOL:
        raise ValueError(f'the pool has {len(ids):,} rows, {len(np.unique(ids)):,} distinct; the paper\'s had {POOL:,} '
                         'distinct sequences')
    keep = ~np.isin(ids, np.array(sorted(excluded), dtype='S64'))
    substrate = np.sort(ids[keep])
    log(f'pool: {len(ids):,} sequences, {int((~keep).sum()):,} of them HQ, development or test, {len(substrate):,} drawn from')
    order = np.random.default_rng(SEED).permutation(len(substrate))
    return [h.decode() for h in substrate[order[:DRAW]]]


def vocabulary_rule(ontology, terms):
    """True if every is_a ancestor of every term (MF root excepted) is in the vocabulary; terms must be live MF ids below
    the root."""
    for t in terms:
        term = ontology.terms.get(t)
        if t == ROOT_TERM or term is None or term['obsolete'] or term['namespace'] != 'molecular_function':
            raise ValueError(f'annotation with {t}, which is not a live molecular-function term below the root')
        if any(a not in ontology.index and a != ROOT_TERM for a in ontology.ancestors(t)):
            return False
    return True


def from_pool(root, ontology, threads):
    """A new draw from external/computational_pool: (protein ids, ranks, accessions, sequences, annotation rows, labels)."""
    pool = root / 'external' / 'computational_pool'
    excluded = {h.encode() for name in ('hq', 'development', 'test') for h in load(root, name, labels=False).seq_ids}
    drawn = draw_order(pool, excluded)
    wanted = pa.array(drawn)
    table = pads.dataset(pool / 'proteins.parquet').to_table(filter=pc.field('seq_id').isin(wanted)).to_pydict()
    record = {h: (s.upper(), acc) for h, s, acc in zip(table['seq_id'], table['sequence'], table['accessions'])}
    bad = [h for h in drawn if h not in record or seq_id(record[h][0]) != h]
    if bad:
        raise ValueError(f'{len(bad)} drawn sequences missing from the pool or not matching their seq_id, e.g. {bad[0]}')
    assertions = pads.dataset(pool / 'assertions.parquet').to_table(filter=pc.field('seq_id').isin(wanted)).to_pydict()
    rows = defaultdict(list)
    for row in zip(*(assertions[k] for k in ('seq_id', 'go_id', 'evidence', 'assigned_by'))):
        rows[row[0]].append(annotation_row(*row))
    missing = [h for h in drawn if not rows[h]]
    if missing:
        raise ValueError(f'{len(missing)} drawn sequences have no annotation, e.g. {missing[0]}')

    exact, hits = screen(root, {h: record[h][0] for h in drawn}, root / 'selection' / 'hq_random_screen', threads)
    outside = {h for h in drawn if not vocabulary_rule(ontology, {r[1] for r in rows[h]})}
    kept = [(k, h) for k, h in enumerate(drawn, 1) if h not in exact | hits | outside][:KEEP]
    if len(kept) < KEEP:
        raise ValueError(f'only {len(kept):,} of {DRAW:,} draws pass the screen and the vocabulary rule')
    picks = [h for _, h in kept]
    considered = drawn[:kept[-1][0]]
    log(f'hq_random: the first {len(considered):,} draws give {KEEP:,} picks; removed by the screen '
        f'{sum(h in exact | hits for h in considered):,}, by the vocabulary rule {sum(h in outside for h in considered):,}')
    pick_rows = [r for h in picks for r in rows[h]]
    labels, dropped = ontology.label_matrix([(h, g) for h, g, _, _ in pick_rows], {h: i for i, h in enumerate(picks)})
    if dropped:
        raise ValueError(f'{len(dropped)} annotations of the picks map to no live term')
    return (picks, [k for k, _ in kept], [record[h][1].split(',') if record[h][1] else [] for h in picks],
            [record[h][0] for h in picks], pick_rows, labels)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--from-pool', action='store_true', help='draw again from external/computational_pool')
    parser.add_argument('--threads', type=int, default=16, help='DIAMOND threads (--from-pool)')
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    ontology = Ontology(root)
    added = from_pool(root, ontology, a.threads) if a.from_pool else from_release(root, ontology)
    if len(added[0]) != KEEP:
        raise ValueError(f'{len(added[0]):,} added proteins, not {KEEP:,}')
    write_added(root, 'hq_random', ontology, *added, 'added', 'order_rank')


if __name__ == '__main__':
    main()
