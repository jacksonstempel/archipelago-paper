"""Gains by distance from HQ (Supplementary Information): the test proteins sorted by their distance to the nearest HQ
protein and split into five groups of 1,000, with the test micro-AUPRC differences between training sets in each group and
over all of them, for the masked, reviewed and experimental-only test labels.

    python -m retrain.scoring.distance_groups [--data DIR] [--device cuda|cpu]

The distance is that of retrain.scoring.distances, from inputs/embeddings/esm2_15b_l48/{test,hq}; ties keep the benchmark
order. Experimental-only labels keep only the proteins that have such a label. A score is computed over the proteins of the
group, with the masked metrics' code (retrain.scoring.masked): the scorer's 10,000 resamples of the test clusters give the
paired 95% intervals, and a resample that draws every cluster once gives the point estimate. A model's value is the mean over
its seeds.

Writes results/distance_groups/: groups.tsv (protein_id, distance, group; group 0 is the nearest) and gains.tsv (labels, group,
proteins, positives, median_distance, method, training_set, minus_training_set, value, minus_value, difference, ci_low,
ci_high; group `all` is every test protein). Needs results/main. Two minutes on one GPU.
"""
import argparse

import numpy as np
import torch

from archipelago.evaluate import cluster_draws
from archipelago.go import Propagation, read_parents
from archipelago.io import read_matrix
from retrain.common import add_argument, data_root, load, log
from retrain.scoring.distances import nearest, unit
from retrain.scoring.masked import replicate_auprc
from retrain.scoring.score import PUBLISHED, bundle, fits
from retrain.scoring.test_labels import annotation_sets

GROUPS = 5
MODELS = [('finch', s) for s in ('hq', 'hq_random', 'archipelago')] + [(m, s) for m in PUBLISHED for s in ('hq', 'archipelago')]
DIFFERENCES = [(m, 'archipelago', 'hq') for m in ('finch', *PUBLISHED)] + \
              [('finch', 'hq_random', 'hq'), ('finch', 'archipelago', 'hq_random')]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--device', choices=['cuda', 'cpu'], help='default: cuda if available')
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    out = root / 'results' / 'distance_groups'
    out.mkdir(parents=True, exist_ok=False)
    device = torch.device(a.device or ('cuda' if torch.cuda.is_available() else 'cpu'))
    torch.backends.cuda.matmul.allow_tf32 = False
    tensor = lambda x: torch.as_tensor(x, device=device)  # noqa: E731

    test, vocabulary, pairs, labels_of = annotation_sets(root)
    embeddings = root / 'inputs' / 'embeddings' / 'esm2_15b_l48'
    hq_ids = load(root, 'hq', labels=False, sequences=False).seq_ids
    distance = nearest(unit(read_matrix(embeddings / 'test', test.seq_ids)[1], device),
                       unit(read_matrix(embeddings / 'hq', hq_ids)[1], device))
    group = np.empty(len(distance), dtype=int)
    group[np.argsort(distance, kind='stable')] = np.arange(len(distance)) * GROUPS // len(distance)
    with open(out / 'groups.tsv', 'w', newline='\n') as handle:
        handle.write('protein_id\tdistance\tgroup\n')
        handle.writelines(f'{p}\t{d:.5f}\t{g}\n' for p, d, g in zip(test.seq_ids, distance, group))

    official, positive, experimental = labels_of(pairs['all']), labels_of(pairs['other']), labels_of(pairs['exp'])
    has_experimental = experimental.any(1)
    labels = {'masked': (tensor(positive), tensor(~(labels_of(pairs['all_goa']) & ~positive))),
              'reviewed': (tensor(official), tensor(np.ones(official.shape, dtype=bool))),
              'experimental': (tensor(experimental), tensor(has_experimental)[:, None].expand(official.shape))}
    subsets = {'all': np.ones(len(group), dtype=bool), **{str(g): group == g for g in range(GROUPS)}}

    propagate = Propagation(read_parents(bundle(root) / 'go-basic.obo', vocabulary))
    cluster, draws = cluster_draws(test.columns['cluster'])
    cluster = tensor(cluster)
    draws = torch.as_tensor(draws.T, dtype=torch.float64, device=device)
    once = torch.ones((draws.shape[0], 1), dtype=torch.float64, device=device)  # every cluster once: the test set itself

    point, replicates = {}, {}
    for method, training_set, seeds in fits(root, 'main'):
        if (method, training_set) not in MODELS:
            continue
        per_seed = {}
        for files in seeds.values():
            s = propagate(torch.as_tensor(np.load(files['test']), dtype=torch.float32, device=device))
            for version, (y, keep) in labels.items():
                for name, chosen in subsets.items():
                    k = keep & tensor(chosen)[:, None]
                    per_seed.setdefault((version, name), []).append((float(replicate_auprc(s, y, k, cluster, once)[0]),
                                                                     replicate_auprc(s, y, k, cluster, draws)))
        for (version, name), values in per_seed.items():
            point[version, name, method, training_set] = float(np.mean([v for v, _ in values]))
            replicates[version, name, method, training_set] = np.mean([r for _, r in values], 0)
        log(method, training_set)

    with open(out / 'gains.tsv', 'w', newline='\n') as handle:
        handle.write('labels\tgroup\tproteins\tpositives\tmedian_distance\tmethod\ttraining_set\tminus_training_set\tvalue\t'
                     'minus_value\tdifference\tci_low\tci_high\n')
        for version, (y, keep) in labels.items():
            for name, chosen in subsets.items():
                proteins = chosen & (has_experimental if version == 'experimental' else True)
                positives = int((y & keep & tensor(chosen)[:, None]).sum())
                for method, a_set, b_set in DIFFERENCES:
                    first, second = point[version, name, method, a_set], point[version, name, method, b_set]
                    low, high = np.percentile(replicates[version, name, method, a_set] - replicates[version, name, method, b_set],
                                              [2.5, 97.5])
                    handle.write(f'{version}\t{name}\t{int(proteins.sum())}\t{positives}\t{float(np.median(distance[proteins]))!r}\t'
                                 f'{method}\t{a_set}\t{b_set}\t{first!r}\t{second!r}\t{first - second!r}\t{float(low)!r}\t{float(high)!r}\n')
    log(f'wrote {out}')


if __name__ == '__main__':
    main()
