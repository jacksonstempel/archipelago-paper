"""Training cost of Finch and the published methods (Supplementary Information): results/training_cost.tsv with, per method
and training set of the main comparison, the number of runs and their summed training time in GPU-hours.

    python -m retrain.scoring.training_cost [--data DIR]

The time is each fit's recorded `train_seconds`: the time spent training, including the trainer's own evaluations on the
development set, but not reading the inputs or predicting. Every fit trains on one GPU. One run of Finch is its two
networks (runs/finch/<set>/seed<k>/<readout>/history.json); one run of DeepGO-SE is its ensemble, for which ten worlds are
trained (runs/deepgose/<set>/seed0/members/world<k>/record.json); TransFew and DPFunc: each fit's record.json.
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
        parts = [run / readout / 'history.json' for readout in config()['weights']]
    elif method == 'deepgose':
        parts = [run / 'members' / f'world{k}' / 'record.json' for k in range(WORLDS)]
    else:
        parts = [run / 'record.json']
    return sum(json.loads(p.read_text())['train_seconds'] for p in parts)


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
