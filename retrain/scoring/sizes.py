"""Sizes of the training sets and the benchmark (Supplementary Table 2): results/dataset_sizes.tsv with, per dataset, its
proteins, propagated positive pairs and the number of terms with at least one positive.

    python -m retrain.scoring.sizes [--data DIR]
"""
import argparse

from retrain.common import add_argument, data_root, label_counts, load, log

DATASETS = ('original_transfew', 'original_deepgose', 'original_dpfunc', 'goa', 'hq', 'hq_random', 'archipelago',
            'development', 'test')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argument(parser)
    root = data_root(parser.parse_args().data)
    lines = ['dataset\tproteins\tpropagated_pairs\tterms_used']
    for name in DATASETS:
        c = label_counts(load(root, name).labels)
        lines.append(f"{name}\t{c['proteins']}\t{c['propagated_pairs']}\t{c['terms_with_positive']}")
    out = root / 'results' / 'dataset_sizes.tsv'
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w', newline='\n') as handle:
        handle.write('\n'.join(lines) + '\n')
    log(f'wrote {out}')


if __name__ == '__main__':
    main()
