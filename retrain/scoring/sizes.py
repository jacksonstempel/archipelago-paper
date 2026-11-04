"""Sizes of the training sets and the benchmark (Supplementary Table 2): results/dataset_sizes.tsv with, per dataset, its
proteins, propagated positive pairs and the number of terms with at least one positive.

    python -m retrain.scoring.sizes [--data DIR] [--labels reviewed|masked]

--labels masked writes results_masked/dataset_sizes.tsv, whose test row counts the masked test labels (the pairs that the reviewed
GOA assertions add are neither positive nor negative, so they are not counted); the other rows are the same.
"""
import argparse

from retrain.common import add_argument, data_root, label_counts, load, log
from retrain.scoring.score import RESULTS, results_dir
from retrain.scoring.test_labels import annotation_sets, version_labels

DATASETS = ('original_transfew', 'original_deepgose', 'original_dpfunc', 'goa', 'hq', 'hq_random', 'archipelago',
            'development', 'test')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--labels', choices=list(RESULTS), default='reviewed', help='the version of the main test labels')
    add_argument(parser)
    a = parser.parse_args(argv)
    root = data_root(a.data)
    lines = ['dataset\tproteins\tpropagated_pairs\tterms_used']
    for name in DATASETS:
        labels = load(root, name, scoring=True).labels  # the benchmark's labels are counted, not used
        if name == 'test' and a.labels == 'masked':
            _, _, pairs, labels_of = annotation_sets(root)
            labels = version_labels(pairs, labels_of, 'masked')[0]
        c = label_counts(labels)
        lines.append(f"{name}\t{c['proteins']}\t{c['propagated_pairs']}\t{c['terms_with_positive']}")
    out = results_dir(root, a.labels) / 'dataset_sizes.tsv'
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w', newline='\n') as handle:
        handle.write('\n'.join(lines) + '\n')
    log(f'wrote {out}')


if __name__ == '__main__':
    main()
