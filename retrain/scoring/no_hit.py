"""Test proteins without training relatives: every model of the main run scored on the test proteins that have no DIAMOND
hit in any training set at progressively looser E-value cutoffs (Supplementary Information).

    python -m retrain.scoring.no_hit [--data DIR] [--labels reviewed|masked] [--search] [--threads 32] [--device cuda|cpu]

1. The subsets: for each cutoff c in 1e-10, 1e-6, 1e-3, 0.01, 0.1, 1, 10, the test proteins without a hit at E <= c. By
   default they are the paper's, from the released table release/benchmark_no_hit (no_hit_up_to >= c): 4,938, 4,781,
   4,307, 3,346, 2,463, 1,489 and 723 proteins. The paper's search ran on every distinct training protein of the main
   comparisons and of one further method whose original training data (374,070 proteins, 316,033 in no other set) is
   not part of this release, so --search, which repeats it on the sets here (goa, hq, hq_random, archipelago and the
   published methods' original data; blastp --ultra-sensitive, E <= 10,000, at most 25 targets per query, the screen's
   fixed database size; results/no_hit/test_hits.tsv), keeps more test proteins: up to 972 instead of 723 at E <= 10. It
   says how its subsets differ from the released ones.
2. Each subset, listed in scoring/no_hit_E<c>.txt, is scored with the scorer's --proteins (the main run's fits):
   results/no_hit/E<c>/, with contrasts.tsv as retrain.scoring.score writes it. Sizes in results/no_hit/subsets.tsv.
3. results/no_hit/differences.tsv: per subset (`none` = all test proteins, from results/main), model and earlier training
   set, the test micro-AUPRC of ArchipelaGO minus that set (seed means) with its paired 95% cluster-bootstrap interval.

--labels masked scores the same subsets on the masked test labels (the subsets stay the same) into results_masked/no_hit; it needs
results_masked/main (retrain.scoring.score main --labels masked) instead of results/main.

Needs results/main (python -m retrain.scoring.score main); --search needs DIAMOND (2.2.4 in the paper) on PATH or in
$DIAMOND.
"""
import argparse
import csv
import json
import subprocess
from pathlib import Path

import numpy as np

from retrain.common import DIAMOND_DBSIZE, add_argument, data_root, diamond, load, log, write_fasta
from retrain.scoring.score import RESULTS, bundle, contrasts, fits, masked_config, results_dir, run_scorer, write_config

TRAINING_SETS = ('goa', 'hq', 'hq_random', 'archipelago', 'original_deepgose', 'original_dpfunc', 'original_transfew')
CUTOFFS = (1e-10, 1e-6, 1e-3, 1e-2, 1e-1, 1.0, 10.0)
BASES = ('original', 'goa', 'hq', 'hq_random')


def search(root, out, threads):
    """Best E-value of each test protein against the union of the training sets: {seq_id: E}."""
    seen, union = set(), []
    for name in TRAINING_SETS:
        ds = load(root, name, labels=False)
        union += [(h, s) for h, s in zip(ds.seq_ids, ds.sequences) if not (h in seen or seen.add(h))]
    log(f'{len(union):,} distinct training proteins')
    write_fasta(out / 'union.fasta', union)
    test = load(root, 'test', labels=False)
    write_fasta(out / 'test.fasta', zip(test.seq_ids, test.sequences))
    binary = diamond()
    subprocess.run([binary, 'makedb', '--in', out / 'union.fasta', '--db', out / 'union', '--threads', str(threads)],
                   check=True)
    subprocess.run([binary, 'blastp', '--db', out / 'union.dmnd', '--query', out / 'test.fasta', '--out',
                    out / 'test_hits.tsv', '--threads', str(threads), '--ultra-sensitive', '--evalue', '10000',
                    '--max-target-seqs', '25', '--dbsize', str(DIAMOND_DBSIZE), '--outfmt', '6', 'qseqid', 'sseqid',
                    'bitscore', 'evalue'], check=True)
    (out / 'union.fasta').unlink()
    best = {}
    with open(out / 'test_hits.tsv') as handle:
        for query, _, _, evalue in csv.reader(handle, delimiter='\t'):
            best[query] = min(best.get(query, float('inf')), float(evalue))
    return best


def released_subsets(root):
    """{cutoff: the test proteins without a hit at E <= cutoff} from the released table (the paper's subsets)."""
    import pyarrow.parquet as pq
    table = pq.read_table(Path(root) / 'release' / 'benchmark_no_hit' / 'test.parquet',
                          columns=['protein_id', 'no_hit_up_to']).to_pydict()
    listed = {v for v in table['no_hit_up_to'] if v is not None}
    if not listed <= set(CUTOFFS):
        raise ValueError(f'benchmark_no_hit lists cutoffs {sorted(listed)}, expected {CUTOFFS}')
    return {c: {h for h, v in zip(table['protein_id'], table['no_hit_up_to']) if v is not None and v >= c}
            for c in CUTOFFS}


def seed_means(results):
    """({(method, training set): test AUPRC seed mean}, {(method, training set): replicate AUPRC of the seed mean})."""
    with open(results / 'metrics.tsv', newline='') as handle:
        value = {(r['method'], r['training_set']): float(r['value']) for r in csv.DictReader(handle, delimiter='\t')
                 if r['seed'] == 'mean' and r['split'] == 'test' and r['metric'] == 'auprc'}
    z = np.load(results / 'replicates.npz')
    column = [str(m) for m in z['metrics']].index('auprc')
    replicates = {tuple(k.split('/')[:2]): z[k][:, column] for k in z.files if k.endswith('/mean')}
    return value, replicates


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--search', action='store_true',
                        help="search the training sets here instead of reading the paper's subsets (see above)")
    parser.add_argument('--threads', type=int, default=32)
    parser.add_argument('--device', choices=['cuda', 'cpu'])
    parser.add_argument('--labels', choices=list(RESULTS), default='reviewed', help='the version of the main test labels')
    add_argument(parser)
    a = parser.parse_args(argv)
    root = data_root(a.data)
    out = results_dir(root, a.labels) / 'no_hit'
    out.mkdir(parents=True, exist_ok=False)
    released = released_subsets(root) if (root / 'release' / 'benchmark_no_hit').exists() else None
    if a.search:
        best = search(root, out, a.threads)
        subsets_of = {c: {h for h in load(root, 'test', labels=False, sequences=False).seq_ids
                          if best.get(h, float('inf')) > c} for c in CUTOFFS}
        for c in CUTOFFS:
            log(f"E <= {c:g}: {len(subsets_of[c]):,} test proteins without a hit here"
                + ('' if released is None else f", {len(released[c]):,} in the paper's subset, "
                   f"{len(subsets_of[c] - released[c]):,} not in it, {len(released[c] - subsets_of[c]):,} of it missing"))
    elif released is None:
        raise SystemExit("release/benchmark_no_hit is missing: download the dataset (data/README.md) or pass --search")
    else:
        subsets_of = released
        log("the paper's subsets, from release/benchmark_no_hit")
    if a.labels == 'masked':
        config = masked_config(root, 'main', 'no_hit')
        data = (config.parent / json.loads(config.read_text())['data']).resolve()
    else:
        data = bundle(root)
        config = write_config(root / 'scoring' / 'no_hit.json', data, fits(root, 'main'))
    with open(data / 'benchmark' / 'test.tsv', newline='') as handle:
        table = list(csv.DictReader(handle, delimiter='\t'))
    subsets = [('none', len(table), len({r['cluster'] for r in table}), results_dir(root, a.labels) / 'main')]
    for cutoff in CUTOFFS:
        keep = {r['protein_id'] for r in table if r['protein_id'] in subsets_of[cutoff]}
        clusters = len({r['cluster'] for r in table if r['protein_id'] in keep})
        name = f'{cutoff:g}'
        log(f'E <= {name}: {len(keep):,} test proteins without a hit, {clusters:,} clusters')
        listed = root / 'scoring' / f'no_hit_E{name}.txt'
        listed.write_text(''.join(f'{h}\n' for h in sorted(keep)), encoding='utf-8')
        run_scorer(root, config, out / f'E{name}', a.device, labels=a.labels, proteins=listed)
        contrasts(out / f'E{name}')
        subsets.append((name, len(keep), clusters, out / f'E{name}'))
    with open(out / 'subsets.tsv', 'w', newline='\n') as handle:
        handle.write('cutoff\tproteins\tclusters\n')
        handle.writelines(f'{float(c)}\t{n}\t{k}\n' for c, n, k, _ in subsets[1:])
    with open(out / 'differences.tsv', 'w', newline='\n') as handle:
        handle.write('subset\tproteins\tmethod\tbase\tdifference\tci_low\tci_high\n')
        for name, n, _, results in subsets:
            value, replicates = seed_means(results)
            for (method, training_set) in sorted(value):
                if training_set != 'archipelago':
                    continue
                for base in BASES:
                    if (method, base) in value:
                        low, high = np.percentile(replicates[method, 'archipelago'] - replicates[method, base], [2.5, 97.5])
                        handle.write(f"{name}\t{n}\t{method}\t{base}\t{value[method, 'archipelago'] - value[method, base]!r}"
                                     f'\t{float(low)!r}\t{float(high)!r}\n')
    log(f'wrote {out / "differences.tsv"}')


if __name__ == '__main__':
    main()
