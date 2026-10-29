"""Test proteins without training relatives: every model of the main run scored on the test proteins that have no DIAMOND
hit in any training set at progressively looser E-value cutoffs (Supplementary Information).

    python -m retrain.scoring.no_hit [--data DIR] [--threads 32] [--device cuda|cpu]

1. One DIAMOND database of every distinct training protein of the main comparisons (goa, hq, hq_random, archipelago and the
   published methods' original data); the 5,000 test proteins searched against it (blastp --ultra-sensitive,
   E <= 10,000, at most 25 targets per query, the screen's fixed database size): results/no_hit/test_hits.tsv.
2. For each cutoff c in 1e-10, 1e-6, 1e-3, 0.01, 0.1, 1, 10: the test proteins whose best hit has E > c (or that have no
   hit), scored with the scorer (the main run's fits; test labels and clusters of those proteins only):
   results/no_hit/E<c>/, with contrasts.tsv as retrain.scoring.score writes it. Sizes in
   results/no_hit/subsets.tsv.
3. results/no_hit/differences.tsv: per subset (`none` = all test proteins, from results/main), model and earlier training
   set, the test micro-AUPRC of ArchipelaGO minus that set (seed means) with its paired 95% cluster-bootstrap interval.

Needs results/main (python -m retrain.scoring.score main) and DIAMOND (2.2.4 in the paper) on PATH or in $DIAMOND.
The paper's union also held the original training data of one further method that is not part of this release, so the
subsets built here can keep a few test proteins more than the paper's.
"""
import argparse
import csv
import subprocess
from pathlib import Path

import numpy as np

from retrain.common import DIAMOND_DBSIZE, add_argument, data_root, diamond, load, log, write_fasta
from retrain.scoring.score import bundle, contrasts, fits, run_scorer, test_subset

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


def seed_means(results):
    """({(method, training set): test AUPRC seed mean}, {(method, training set): replicate AUPRC of the seed mean})."""
    with open(results / 'metrics.tsv', newline='') as handle:
        value = {(r['method'], r['training_set']): float(r['value']) for r in csv.DictReader(handle, delimiter='\t')
                 if r['seed'] == 'mean' and r['split'] == 'test' and r['metric'] == 'auprc'}
    z = np.load(results / 'replicates.npz')
    column = [str(m) for m in z['metrics']].index('auprc')
    replicates = {tuple(k.split('/')[:2]): z[k][:, column] for k in z.files if k.endswith('/mean')}
    return value, replicates


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--threads', type=int, default=32)
    parser.add_argument('--device', choices=['cuda', 'cpu'])
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    out = root / 'results' / 'no_hit'
    out.mkdir(parents=True, exist_ok=False)
    best = search(root, out, a.threads)
    data = bundle(root)
    models = fits(root, 'main')
    with open(data / 'benchmark' / 'test.tsv', newline='') as handle:
        table = list(csv.DictReader(handle, delimiter='\t'))
    with open(data / 'benchmark' / 'test_labels.tsv', newline='') as handle:
        labels = [(r['protein_id'], r['go_id']) for r in csv.DictReader(handle, delimiter='\t')]
    subsets = [('none', len(table), len({r['cluster'] for r in table}), root / 'results' / 'main')]
    for cutoff in CUTOFFS:
        keep = {r['protein_id'] for r in table if best.get(r['protein_id'], float('inf')) > cutoff}
        clusters = len({r['cluster'] for r in table if r['protein_id'] in keep})
        name = f'{cutoff:g}'
        log(f'E <= {name}: {len(keep):,} test proteins without a hit, {clusters:,} clusters')
        config = test_subset(root, data, root / 'scoring' / f'no_hit_E{name}', keep, labels, models)
        run_scorer(root, config, out / f'E{name}', a.device)
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
