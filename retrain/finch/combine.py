"""Finch's scores: per training set and seed, the weighted mean of its readouts' raw scores (weights in configs/finch.json).

    python -m retrain.finch.combine --training-set archipelago [--seeds 9301 9302 9303] [--model finch] [--data DIR]

Reads runs/finch/<training set>/seed<k>/<readout>/{dev,test}_scores.npy (retrain.finch.train) and writes
runs/<model>/<training set>/seed<k>/{dev_scores.npy, test_scores.npy, record.json}. --model finch (the default) is Finch
(ESM2-15B 0.725, ProtT5 0.275); the models of the recipe ablation are listed under 'ladder' in configs/finch.json
(--model ladder_linear_esm2_3b, ...), and Finch's variants under 'variants' (--model finch_prott5_l24). The paper's
expression: the sum, in the listed order of the readouts, of float32 weight x raw float32 score. Every readout must be
complete and be the readout of this training set and seed (its seed plus its seed_offset). Scores stay raw: the scorer
propagates them up the GO graph and forms the seed ensemble (the mean of the seeds' propagated scores).
"""
import argparse
import json

import numpy as np

from retrain.common import add_argument, data_root, fit_folder, log, write_fit
from retrain.finch import config


def combine(root, training_set, seed, model, weights, settings):
    """Write the combination `model` (readout weights) of one training set and seed."""
    readouts = fit_folder(root, 'finch', training_set, seed)
    paths, histories = {}, {}
    for readout in weights:
        path = readouts / readout
        if not all((path / name).is_file() for name in ('readout.pt', 'dev_scores.npy', 'test_scores.npy', 'history.json')):
            raise FileNotFoundError(f'{path}: readout not trained (retrain.finch.train)')
        history = json.loads((path / 'history.json').read_text(encoding='utf-8'))
        expected = (readout, training_set, seed + settings['readouts'][readout].get('seed_offset', 0))
        found = tuple(history.get(k) for k in ('readout', 'training_set', 'seed'))
        if found != expected:
            raise ValueError(f'{path}: history.json is readout, training set, seed {found}, expected {expected}')
        paths[readout], histories[readout] = path, history
    scores = {}
    for split in ('dev', 'test'):
        parts = {m: np.load(paths[m] / f'{split}_scores.npy') for m in weights}
        scores[split] = sum(np.float32(weights[m]) * parts[m] for m in weights).astype(np.float32)
    folder = fit_folder(root, model, training_set, seed)
    write_fit(folder, scores['dev'], scores['test'], dict(
        method=model, training_set=training_set, seed=seed, weights=weights, propagated=False,
        best_step={m: h['best_step'] for m, h in histories.items()},
        members={m: dict(path=str(paths[m]), seed=h['seed'], best_step=h['best_step'], dev_metrics=h['dev_metrics'])
                 for m, h in histories.items()}))
    log(f'wrote {folder}')


def main():
    settings = config()
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--training-set', required=True)
    parser.add_argument('--seeds', nargs='+', type=int, default=settings['seeds'])
    models = {'finch': settings['weights'], **settings['ladder'], **settings['variants']}
    parser.add_argument('--model', choices=list(models), default='finch')
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    for seed in a.seeds:
        combine(root, a.training_set, seed, a.model, models[a.model], settings)


if __name__ == '__main__':
    main()
