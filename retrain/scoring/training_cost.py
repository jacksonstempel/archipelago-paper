"""Training cost of Finch and the published methods (Supplementary Information): results/training_cost.tsv with, per method
and training set of the main comparison, the number of runs and their summed training time in GPU-hours.

    python -m retrain.scoring.training_cost [--data DIR]

The time is what each trainer records, as the paper's table took it; every fit trains on one GPU, and every window
includes the trainer's own evaluations on the development set:
- Finch: one run is its two networks, each from the start of training to the predicted test set (saving and the reload
  check included, reading the inputs not): `train_seconds` of runs/finch/<set>/seed<k>/<readout>/history.json;
- DeepGO-SE: one run is its ensemble, for which ten worlds are trained, each from reading its inputs to the predicted test
  set: the sum of `wall_seconds` of runs/deepgose/<set>/seed0/members/world<k>/record.json;
- TransFew: from reading the inputs to the predicted test set: `wall_seconds` of each fit's record.json;
- DPFunc: from reading the inputs to the end of the last epoch (prediction not included): `train_seconds` of each fit's
  record.json.
The paper's ArchipelaGO Finch runs of one seed shared their GPUs with a second fit, so their recorded time is longer
than a run on its own GPU.
"""
import argparse
import json

from retrain.common import add_argument, data_root, log
from retrain.deepgose.train import WORLDS
from retrain.finch import config
from retrain.scoring.score import PUBLISHED, RUNS


def run_seconds(run, method):
    """The training seconds of one run (a fit folder seed<k>)."""
    if method == 'finch':
        parts, key = [run / readout / 'history.json' for readout in config()['weights']], 'train_seconds'
    elif method == 'deepgose':
        parts, key = [run / 'members' / f'world{k}' / 'record.json' for k in range(WORLDS)], 'wall_seconds'
    else:
        parts, key = [run / 'record.json'], {'transfew': 'wall_seconds', 'dpfunc': 'train_seconds'}[method]
    return sum(json.loads(p.read_text())[key] for p in parts)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argument(parser)
    root = data_root(parser.parse_args().data)
    lines = ['method\ttraining_set\truns\tgpu_hours']
    for method in ('finch', *PUBLISHED):
        for training_set in RUNS['main'][method]:
            runs = sorted((root / 'runs' / method / training_set).glob('seed*'))
            if not runs:
                raise FileNotFoundError(f'no fit of {method} on {training_set}')
            hours = sum(run_seconds(run, method) for run in runs) / 3600
            lines.append(f'{method}\t{training_set}\t{len(runs)}\t{hours!r}')
            log(f'{method} on {training_set}: {len(runs)} runs, {hours:.2f} GPU-hours')
    out = root / 'results' / 'training_cost.tsv'
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w', newline='\n') as handle:
        handle.write('\n'.join(lines) + '\n')
    log(f'wrote {out}')


if __name__ == '__main__':
    main()
