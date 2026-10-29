"""Finch's scores: per training set and seed, the weighted mean of its readouts' raw scores (weights in configs/finch.json).

    python -m retrain.finch.combine --training-set archipelago [--seeds 9301 9302 9303] [--model finch] [--data DIR]

Reads runs/finch/<training set>/seed<k>/<readout>/{dev,test}_scores.npy (retrain.finch.train) and writes
runs/<model>/<training set>/seed<k>/{dev_scores.npy, test_scores.npy, record.json}. --model finch (the default) is Finch
(ESM2-15B 0.725, ProtT5 0.275); the models of the recipe ablation are listed under 'ladder' in configs/finch.json
(--model ladder_linear_esm2_3b, ...). Scores stay raw: the scorer propagates them up the GO graph and forms the seed
ensemble (the mean of the seeds' propagated scores).
"""
import argparse
import json

import numpy as np

from archipelago.finch import combine
from retrain.common import add_argument, data_root, fit_folder, log, write_fit
from retrain.finch import config


def main():
    settings = config()
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--training-set', required=True)
    parser.add_argument('--seeds', nargs='+', type=int, default=settings['seeds'])
    parser.add_argument('--model', choices=['finch', *settings['ladder']], default='finch')
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    weights = settings['weights'] if a.model == 'finch' else settings['ladder'][a.model]
    for seed in a.seeds:
        readouts = fit_folder(root, 'finch', a.training_set, seed)
        scores, members = {}, {}
        for readout in weights:
            member = readouts / readout
            if not (member / 'readout.pt').is_file() or not (member / 'test_scores.npy').is_file():
                raise FileNotFoundError(f'{member}: readout not trained (retrain.finch.train)')
            members[readout] = json.loads((member / 'history.json').read_text(encoding='utf-8'))['best_step']
            scores[readout] = {split: np.load(member / f'{split}_scores.npy') for split in ('dev', 'test')}
        dev = combine({r: s['dev'] for r, s in scores.items()}, weights)
        test = combine({r: s['test'] for r, s in scores.items()}, weights)
        folder = fit_folder(root, a.model, a.training_set, seed)
        write_fit(folder, dev, test, dict(method=a.model, training_set=a.training_set, seed=seed, weights=weights,
                                          best_step=members, propagated=False))
        log(f'wrote {folder}')


if __name__ == '__main__':
    main()
