"""Naive baseline: every protein gets, for every term, the term's frequency among the training proteins.

    python -m retrain.baselines.naive --training-sets goa hq archipelago [--data DIR]

Seconds per training set.
"""
import argparse

import numpy as np

from retrain.baselines import Sets
from retrain.common import add_argument, data_root, log


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--training-sets', nargs='+', required=True)
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    for name in a.training_sets:
        sets = Sets(root, name)
        frequency = np.asarray(sets.train.labels.mean(axis=0), dtype=np.float64).ravel().astype(np.float32)
        sets.write('naive', np.tile(frequency, (len(sets.dev.seq_ids), 1)), np.tile(frequency, (len(sets.test.seq_ids), 1)))
        log(f'naive on {name}: written')


if __name__ == '__main__':
    main()
