"""Train Finch's readouts on one training set and seed, and predict the development and test sets.

    python -m retrain.finch.train --training-set archipelago --seed 9301 [--readouts esm2_15b_l48 prott5_l21]
        [--data DIR] [--device cuda]

Inputs under the data root: datasets/<training set>, datasets/development (labels for early stopping), datasets/test
(proteins only), and the embeddings inputs/embeddings/<readout>/{<training set>,development,test}.{npy,ids.txt}
(scripts/embed.py). The recipe of each readout is configs/finch.json (retrain.finch.recipe; the paper's runs used a
fused-kernel implementation of the same recipe, so retrained scores match up to seed and precision noise, not bit for
bit).

For each readout, writes runs/finch/<training set>/seed<k>/<readout>/: readout.pt (state_dict, readout, seed, recipe,
input_dim, terms, vocabulary, best_step: the format scripts/export_package.py reads), dev_scores.npy and test_scores.npy
(raw sigmoid scores, rows in development/test order) and history.json. The test set is predicted only after the
reloaded checkpoint reproduces the development scores kept during training. One readout takes about 5 to 30 minutes on
one H100 (ESM2-15B on ArchipelaGO is the longest); training proteins and labels are held on the GPU. Host memory: the
dense uint8 label matrix and the FP16 ESM2-15B embeddings, each held twice while copied to the GPU, peak about 26 GB for
archipelago_600k (687,797 proteins) and 9 GB for ArchipelaGO.
"""
import argparse
import json
import time

import numpy as np
import torch

from archipelago.io import read_matrix
from archipelago.readout import Readout, predict
from retrain.finch.recipe import train
from retrain.common import Ontology, add_argument, check_reproduced, data_root, fit_folder, load, log
from retrain.finch import config


def recipe_of(settings, readout):
    recipe = {**settings['training'], **settings['readouts'][readout]}
    recipe.pop('encoder')
    return recipe


def fit(root, training_set, seed, readout, recipe, device):
    ontology = Ontology(root)
    data = load(root, training_set)
    dev = load(root, 'development')
    test = load(root, 'test', labels=False)
    overlap = set(data.seq_ids) & (set(dev.seq_ids) | set(test.seq_ids))
    if overlap:
        raise ValueError(f'{training_set}: {len(overlap)} training proteins are development or test proteins')
    for d in (data, dev):
        if d.labels.shape[1] != len(ontology.vocabulary):
            raise ValueError(f'{d.name}: {d.labels.shape[1]} label columns, the vocabulary has {len(ontology.vocabulary)}')
    stem = root / 'inputs' / 'embeddings' / readout
    _, x = read_matrix(stem / training_set, data.seq_ids)
    _, dev_x = read_matrix(stem / 'development', dev.seq_ids)
    _, test_x = read_matrix(stem / 'test', test.seq_ids)
    for name, m in ((training_set, x), ('development', dev_x), ('test', test_x)):
        if not all(np.isfinite(m[a:a + 65536]).all() for a in range(0, len(m), 65536)):
            raise ValueError(f'{readout}: non-finite {name} embeddings')
    log(f'{readout} on {training_set}, seed {seed}: {len(data.seq_ids):,} training proteins, '
        f'{x.shape[1]} dimensions, {data.labels.shape[1]} terms')
    started = time.monotonic()
    result = train(x, data.labels.toarray(), dev_x, dev.labels.toarray(), recipe, seed, device=device, log=log)

    out = fit_folder(root, 'finch', training_set, seed) / readout
    out.mkdir(parents=True, exist_ok=True)
    torch.save(dict(state_dict=result['state'], readout=readout, seed=seed, recipe=recipe, input_dim=int(x.shape[1]),
                    terms=len(ontology.vocabulary), vocabulary=ontology.vocabulary, best_step=result['best_step']),
               out / 'readout.pt')
    saved = torch.load(out / 'readout.pt', map_location='cpu', weights_only=True)
    model = Readout(saved['input_dim'], recipe['width'], saved['terms'], recipe['dropout'])
    model.load_state_dict(saved['state_dict'])
    model.to(device)
    difference = check_reproduced(result['dev_scores'], predict(model, dev_x, device=device))
    np.save(out / 'dev_scores.npy', result['dev_scores'].astype(np.float32))
    np.save(out / 'test_scores.npy', predict(model, test_x, device=device).astype(np.float32))
    (out / 'history.json').write_text(json.dumps(dict(
        best_step=result['best_step'], steps=result['steps'], reproduced_dev_max_abs=difference,
        seconds=round(time.monotonic() - started), checks=result['history']), indent=1) + '\n', encoding='utf-8')
    log(f"wrote {out}: best development check at update {result['best_step']:,} of {result['steps']:,}")


def main():
    settings = config()
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--training-set', required=True, help='a dataset name under datasets/')
    parser.add_argument('--seed', type=int, required=True, help=f"the paper's seeds: {settings['seeds']}")
    parser.add_argument('--readouts', nargs='+', choices=list(settings['readouts']), default=list(settings['readouts']))
    parser.add_argument('--device', default='cuda' if torch.cuda.is_available() else 'cpu')
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    for readout in a.readouts:
        fit(root, a.training_set, a.seed, readout, recipe_of(settings, readout), a.device)


if __name__ == '__main__':
    main()
