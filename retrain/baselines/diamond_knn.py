"""DIAMOND-kNN baseline (DeepGOPlus's DIAMOND component): development and test proteins are searched against the training
proteins (DIAMOND blastp, at most 25 targets per query, DIAMOND's default that DeepGOPlus keeps); a term's score is the
summed bitscore of the hits carrying it divided by the summed bitscore of all hits. A protein without a hit scores 0.

    python -m retrain.baselines.diamond_knn --training-sets goa hq archipelago [--threads 16] [--data DIR]
    python -m retrain.baselines.diamond_knn --training-sets goa hq archipelago --tune

Two fits per training set:
  diamond_knn         the paper's DIAMOND-kNN (Table 1): --ultra-sensitive search up to E = 1e4, hits kept up to the
                      E-value cutoff (--cutoff, default 1e4: the paper's choice, i.e. no cutoff in effect). The paper chose
                      it on the development set from 1e-6, 1e-5, ..., 1e4 as the best mean development micro-AUPRC over
                      goa, hq and archipelago; --tune repeats that choice (only with exactly those three training sets) and
                      writes the grid to runs/diamond_knn/grid.tsv.
  diamond_knn_screen  --very-sensitive search at the homology screen's own cutoff, E <= 1e-3, which finds almost nothing
                      because the screen removed exactly those hits (the paper's Methods quote its 0.008 micro-AUPRC).
Hits and DIAMOND logs are kept in runs/<method>/<training set>/seed0/work/. Minutes per training set with 16 threads.
"""
import argparse
import csv
import subprocess

import numpy as np
import scipy.sparse as sp

from retrain.baselines import SEED, Sets, share
from retrain.common import add_argument, data_root, diamond, fit_folder, log

GRID = tuple(float(f'1e{k}') for k in range(-6, 5))  # exact decimal values, so --cutoff 1e4 matches
TUNED_OVER = ('goa', 'hq', 'archipelago')
SEARCHES = {'diamond_knn': ['--ultra-sensitive', '--evalue', str(max(GRID))],
            'diamond_knn_screen': ['--very-sensitive', '--evalue', '1e-3']}


def search(sets, method, threads):
    """{split: (rows, columns, bitscores, E-values, shape)}: one hit per query-target pair, the best-scoring one."""
    work = fit_folder(sets.root, method, sets.name, SEED) / 'work'
    work.mkdir(parents=True, exist_ok=True)
    binary, out = diamond(), {}
    with open(work / 'diamond.log', 'w', encoding='utf-8') as log_file:
        subprocess.run([binary, 'makedb', '--in', sets.train.folder / 'sequences.fasta', '--db', work / 'train',
                        '--threads', str(threads)], check=True, stdout=log_file, stderr=subprocess.STDOUT)
        targets = {h: i for i, h in enumerate(sets.train.seq_ids)}
        for split, data in (('dev', sets.dev), ('test', sets.test)):
            result = work / f'{split}_hits.tsv'
            subprocess.run([binary, 'blastp', '--db', work / 'train.dmnd', '--query', data.folder / 'sequences.fasta',
                            '--out', result, '--threads', str(threads), *SEARCHES[method], '--max-target-seqs', '25',
                            '--outfmt', '6', 'qseqid', 'sseqid', 'bitscore', 'evalue'],
                           check=True, stdout=log_file, stderr=subprocess.STDOUT)
            queries = {h: i for i, h in enumerate(data.seq_ids)}
            best = {}
            with open(result, encoding='utf-8') as handle:
                for query, subject, bits, evalue in csv.reader(handle, delimiter='\t'):
                    key = (queries[query], targets[subject])
                    if key not in best or float(bits) > best[key][0]:
                        best[key] = (float(bits), float(evalue))
            keys = np.array(list(best), dtype=np.int64).reshape(-1, 2)
            values = np.array(list(best.values()), dtype=np.float64).reshape(-1, 2)
            out[split] = (keys[:, 0], keys[:, 1], values[:, 0], values[:, 1], (len(data.seq_ids), len(targets)))
            log(f'{method} {sets.name} {split}: {len(best):,} hits for {len(set(keys[:, 0].tolist())):,} proteins')
    return out


def scores(hits, labels, cutoff):
    rows, columns, bits, evalues, shape = hits
    keep = evalues <= cutoff
    return share(sp.csr_matrix((bits[keep], (rows[keep], columns[keep])), shape=shape), labels)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--training-sets', nargs='+', required=True)
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--cutoff', type=float, default=max(GRID), choices=GRID, help='E-value cutoff (default: the paper\'s)')
    parser.add_argument('--tune', action='store_true', help='choose the cutoff on development, as the paper did')
    add_argument(parser)
    a = parser.parse_args()
    if a.tune and sorted(a.training_sets) != sorted(TUNED_OVER):
        raise SystemExit(f'--tune chooses one cutoff over the paper\'s training sets {TUNED_OVER}; pass exactly those')
    root = data_root(a.data)
    searched, grid = {}, []
    for name in a.training_sets:
        sets = Sets(root, name)
        hits = search(sets, 'diamond_knn_screen', a.threads)
        sets.write('diamond_knn_screen', scores(hits['dev'], sets.train.labels, 1e-3),
                   scores(hits['test'], sets.train.labels, 1e-3), search=SEARCHES['diamond_knn_screen'])
        hits = search(sets, 'diamond_knn', a.threads)
        if a.tune:
            for cutoff in GRID:
                grid.append(dict(training_set=name, evalue=cutoff,
                                 dev_auprc=sets.dev_auprc(scores(hits['dev'], sets.train.labels, cutoff))))
                log(f"grid {name} E <= {cutoff:g}: development micro-AUPRC {grid[-1]['dev_auprc']:.4f}")
        searched[name] = (sets, hits)
    cutoff, rule = a.cutoff, 'given (--cutoff)'
    if a.tune:
        with open(root / 'runs' / 'diamond_knn' / 'grid.tsv', 'w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(grid[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(grid)
        auprc = {(g['training_set'], g['evalue']): g['dev_auprc'] for g in grid}
        cutoff = max(GRID, key=lambda c: np.mean([auprc[name, c] for name in searched]))
        rule = f'the best mean development micro-AUPRC over {", ".join(TUNED_OVER)}'
    log(f'E-value cutoff {cutoff:g}: {rule}')
    for name, (sets, hits) in searched.items():
        sets.write('diamond_knn', scores(hits['dev'], sets.train.labels, cutoff),
                   scores(hits['test'], sets.train.labels, cutoff), search=SEARCHES['diamond_knn'], evalue_cutoff=cutoff,
                   cutoff_rule=rule)


if __name__ == '__main__':
    main()
