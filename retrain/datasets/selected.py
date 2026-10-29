"""Training sets that add computational candidates (datasets/candidates) to HQ, each by one selection rule.

    python -m retrain.datasets.selected check-archipelago [--data DIR]
    python -m retrain.datasets.selected dose [--data DIR]
    python -m retrain.datasets.selected controls [--data DIR] [--threads N]

check-archipelago  ArchipelaGO is HQ plus the first 150,000 proteins of the coverage ranking: checks that the ranking and
                   datasets/candidates reproduce the released datasets/archipelago (the same proteins in the same order,
                   with the same labels). Writes nothing.
dose               archipelago_300k, archipelago_600k: HQ plus the first 300,000 or 600,000 proteins of the coverage ranking
                   (rows up to 237,797 are ArchipelaGO's, checked). As for ArchipelaGO, the homology screen is the one the
                   candidates already passed.
controls           The sets of Fig. 3b that add the ranking's ingredients one at a time, 150,000 proteins each:
                     pool_random          a uniformly random order of all candidates: numpy default_rng(20261015)
                                          .permutation over the candidates' rows, which must be in the released pool
                                          order (8,952,628 candidates, checked)
                     filter_random        the same order restricted to the candidates that pass the rare-function filter
                                          (a label held by fewer than 1,000 HQ proteins; 5,287,663 candidates, checked)
                     coverage_unfiltered  the coverage ranking computed without the filter
                   The first 165,000 candidates of each order are screened again against the benchmark (exact sequence or
                   DIAMOND hit, retrain.common.screen); a candidate with a hit is skipped and its order continues.
                   HQ + random computational proteins (hq_random) is built by retrain.datasets.random_draw.

Inputs: datasets/hq, development, test and candidates (retrain.datasets.candidates); the coverage rankings
coverage_ranking.tsv (600,000 rows) and coverage_ranking_unfiltered.tsv (500,000 rows), the paper's own as release files
under release/selection/ (used when present; their lengths are checked) or recomputed by scripts/rank_coverage.py under
selection/ (with a warning: tied candidates can come in another order, so sets beyond ArchipelaGO can differ).
"""
import argparse

import numpy as np

from retrain.common import Ontology, add_argument, data_root, load, log, screen
from retrain.datasets.added import Candidates, ranking, rare_function, sequences, write_added

RANDOM_SEED = 20261015
KEEP, DRAW = 150_000, 165_000
DOSES = {'archipelago_300k': 300_000, 'archipelago_600k': 600_000}
POOL, ELIGIBLE = 8_952_628, 5_287_663  # the paper's candidates and those passing the rare-function filter


def check_archipelago(root, candidates, order):
    archipelago = load(root, 'archipelago')
    n_hq = archipelago.columns['part'].count('hq')
    selected = archipelago.seq_ids[n_hq:]
    if order[:len(selected)] != selected:
        raise ValueError('the coverage ranking does not begin with the released ArchipelaGO selected proteins')
    differ = (candidates.labels[candidates.rows(selected)].astype(np.int8)
              - archipelago.labels[n_hq:].astype(np.int8)).count_nonzero()
    if differ:
        raise ValueError(f'candidate labels of the selected proteins differ from ArchipelaGO in {differ} pairs')
    log(f'ArchipelaGO = HQ + the first {len(selected):,} proteins of the coverage ranking, labels equal: ok')
    return n_hq + len(selected)


def dose(root, ontology, candidates):
    order = ranking(root, 'coverage_ranking.tsv')
    rows = check_archipelago(root, candidates, order)
    for name, n in DOSES.items():
        if len(order) < n:
            raise ValueError(f'the coverage ranking has only {len(order):,} rows; {name} needs {n:,}')
        ids = order[:n]
        accessions, seqs, annotation_rows = candidates.details(ids)
        write_added(root, name, ontology, ids, range(1, n + 1), accessions, seqs, annotation_rows,
                    candidates.labels[candidates.rows(ids)], 'selected', 'selection_rank')
        if load(root, name, labels=False).seq_ids[:rows] != load(root, 'archipelago', labels=False).seq_ids:
            raise ValueError(f'{name}: the first {rows:,} rows are not ArchipelaGO')


def controls(root, ontology, candidates, threads):
    if len(candidates.ids) != POOL:
        raise ValueError(f'{len(candidates.ids):,} candidates; the paper\'s random orders permute its {POOL:,}')
    hq_counts = np.asarray(load(root, 'hq').labels.sum(axis=0)).ravel().astype(np.int64)
    eligible = rare_function(candidates.labels, hq_counts)
    if int(eligible.sum()) != ELIGIBLE:
        raise ValueError(f'{int(eligible.sum()):,} candidates pass the rare-function filter; the paper had {ELIGIBLE:,}')
    permutation = np.random.default_rng(RANDOM_SEED).permutation(len(candidates.ids))
    orders = {'pool_random': [candidates.ids[i] for i in permutation[:DRAW]],
              'filter_random': [candidates.ids[i] for i in permutation[eligible[permutation]][:DRAW]],
              'coverage_unfiltered': ranking(root, 'coverage_ranking_unfiltered.tsv')[:DRAW]}
    drawn = {h for order in orders.values() for h in order}
    exact, hits = screen(root, sequences(candidates.folder / 'sequences.fasta', drawn),
                         root / 'selection' / 'controls_screen', threads)
    related = exact | hits
    for name, order in orders.items():
        kept = [(k, h) for k, h in enumerate(order, 1) if h not in related][:KEEP]
        if len(kept) < KEEP:
            raise ValueError(f'{name}: only {len(kept):,} of the first {len(order):,} candidates pass the screen')
        log(f'{name}: {kept[-1][0] - KEEP} of the first {kept[-1][0]:,} candidates skipped by the screen')
        ids = [h for _, h in kept]
        accessions, seqs, annotation_rows = candidates.details(ids)
        write_added(root, name, ontology, ids, [k for k, _ in kept], accessions, seqs, annotation_rows,
                    candidates.labels[candidates.rows(ids)], 'added', 'order_rank')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('stage', choices=['check-archipelago', 'dose', 'controls'])
    parser.add_argument('--threads', type=int, default=16, help='DIAMOND threads (controls)')
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    candidates = Candidates(root)
    if a.stage == 'check-archipelago':
        check_archipelago(root, candidates, ranking(root, 'coverage_ranking.tsv'))
    elif a.stage == 'dose':
        dose(root, Ontology(root), candidates)
    else:
        controls(root, Ontology(root), candidates, a.threads)


if __name__ == '__main__':
    main()
