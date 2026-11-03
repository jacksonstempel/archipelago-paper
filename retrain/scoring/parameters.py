"""Trainable parameters of the models trained on ArchipelaGO (the SI tables of the recipe ablation and further models):
results/parameters.json, {model: {part: parameters}}; a model's count is the sum of its parts.

    python -m retrain.scoring.parameters [--data DIR]      (CPU; needs the deepgose, transfew and dpfunc extras)

Counted from the checkpoint of each model's first ArchipelaGO run; the frozen protein language models are not counted.
- Finch, the recipe ablation's models and Finch's variant (configs/finch.json): one part per readout, the entries of
  runs/finch/archipelago/seed<first seed>/<readout>/readout.pt (all trained), checked against the formula of its shape.
- The published methods are rebuilt from their model code with the checkpoint loaded strictly; buffers (BatchNorm
  statistics, index tensors) are not parameters.
- DeepGO-SE: one part per world of its ensemble (protein MLP and the GO-Plus class and relation embeddings, all
  trained).
- TransFew (last epoch) and DPFunc (last epoch; its prediction averages the last three epochs of the same network): one
  part per top-level module. TransFew's label GCN (label_embedding), pretrained on GO alone and frozen in training, is
  counted; the sentence encoder that embedded the term definitions for it is not part of the model.
"""
import argparse
import json

import torch

from retrain.common import add_argument, data_root, fit_folder, log
from retrain.finch import config


def modules(model):
    """{top-level module: parameters} of a model."""
    out = {}
    for name, p in model.named_parameters():
        out[name.split('.')[0]] = out.get(name.split('.')[0], 0) + p.numel()
    return out


def readout(root, name, seed):
    """Every entry of a saved readout (LayerNorm, linear layers) is a trained parameter; checked against the formula
    2d + d*h + h + h*T + T (hidden layer of width h) or 2d + d*T + T (linear readout)."""
    path = fit_folder(root, 'finch', 'archipelago', seed) / name / 'readout.pt'
    saved = torch.load(path, map_location='cpu', weights_only=True)
    state, dim, terms = saved['state_dict'], saved['input_dim'], saved['terms']
    hidden = 'network.4.weight' in state
    keys = {'network.0.weight', 'network.0.bias', 'network.1.weight', 'network.1.bias'} | (
        {'network.4.weight', 'network.4.bias'} if hidden else set())
    if set(state) != keys:
        raise ValueError(f'{path}: unexpected readout {sorted(state)}')
    width = state['network.1.weight'].shape[0] if hidden else None
    formula = 2 * dim + (dim * width + width + width * terms + terms if hidden else dim * terms + terms)
    count = sum(v.numel() for v in state.values())
    if count != formula or state['network.0.weight'].shape[0] != dim:
        raise ValueError(f'{path}: {count} parameters, formula {formula}')
    return count


def deepgose(root, terms):
    from retrain.deepgose.model import DeepGOModel
    from retrain.deepgose.train import world_folder
    record = json.loads((fit_folder(root, 'deepgose', 'archipelago', 0) / 'record.json').read_text())
    out = {}
    for k in record['ensemble_worlds']:
        state = torch.load(world_folder(root, 'archipelago', k) / 'best.pt', map_location='cpu')
        model = DeepGOModel(state['net.0.linear.weight'].shape[1], terms, state['go_embed.weight'].shape[0] - terms,
                            state['rel_embed.weight'].shape[0] - 1, 'cpu')
        model.load_state_dict(state)
        out[f'world{k}'] = sum(modules(model).values())
    return out


def transfew(root):
    from retrain.transfew.train import build
    state = torch.load(fit_folder(root, 'transfew', 'archipelago', 0) / 'checkpoint_last.pt', map_location='cpu')
    freq, rare = state['freq_indicies'], state['rare_indicies']
    model = build(root, len(freq) + len(rare), freq, rare)
    model.load_state_dict(state)
    return modules(model)


def dpfunc(root, terms):
    from retrain.dpfunc.train import EPOCHS, new_model
    model = new_model(terms, 'cpu')
    model.load_state_dict(torch.load(fit_folder(root, 'dpfunc', 'archipelago', 0) / f'epoch{EPOCHS}.pt',
                                     map_location='cpu')['model_state_dict'])
    return modules(model)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argument(parser)
    root = data_root(parser.parse_args().data)
    settings = config()
    terms = len(json.loads((root / 'datasets' / 'ontology' / 'vocabulary.json').read_text(encoding='utf-8')))
    models = {'finch': settings['weights'], **settings['ladder'], **settings['variants']}
    counts = {name: readout(root, name, settings['seeds'][0]) for name in {r for m in models.values() for r in m}}
    out = {model: {r: counts[r] for r in readouts} for model, readouts in models.items()}
    out.update(deepgose=deepgose(root, terms), transfew=transfew(root), dpfunc=dpfunc(root, terms))
    path = root / 'results' / 'parameters.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1) + '\n', encoding='utf-8')
    for model, parts in out.items():
        log(f'{model}: {sum(parts.values()) / 1e6:,.1f} million')
    log(f'wrote {path}')


if __name__ == '__main__':
    main()
