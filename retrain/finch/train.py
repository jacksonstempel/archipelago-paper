"""Train Finch's readouts on one training set and seed, and predict the development and test sets.

    python -m retrain.finch.train --training-set archipelago --seed 9301 [--readouts esm2_15b_l48 prott5_l21]
        [--data DIR] [--allow-shared-proteins]

The readouts and their recipes are configs/finch.json: by default Finch's two; the recipe ablation's readouts (Fig. 4a)
and its variant's are trained on ArchipelaGO with --readouts esm2_3b_l36_linear esm2_15b_l48_linear esm2_15b_l48_narrow
esm2_15b_l48_adamw esm2_15b_l48_second esm2_3b_l36_adamw esm1b_l31_adamw prott5_l24. Training is the paper's engine
(retrain/finch/engine.py) and the steps around it are the paper's fit: training proteins and labels are held on the GPU,
the best raw development micro-AUPRC checkpoint is kept, the reloaded checkpoint must reproduce its development scores
(within 1e-5) before the test set is predicted. Needs one CUDA GPU (an H100 80 GB in the paper) and Triton.

Inputs under the data root: datasets/<training set>, datasets/development (labels for checkpoint selection),
datasets/test (proteins only), datasets/ontology, and each readout's FP16 embeddings
inputs/embeddings/<embedding>/{<training set>,development,test}.{npy,ids.txt} (scripts/embed.py, or the release's).

For each readout, writes runs/finch/<training set>/seed<k>/<readout>/: readout.pt (state_dict, readout, seed, recipe,
input_dim, terms, vocabulary, best_step: the format scripts/export_package.py reads), dev_scores.npy and test_scores.npy
(raw sigmoid scores, float32, rows in development/test order) and history.json (updates, the selected update, the
development metrics of the checkpoint, the reload check, train_seconds and every development check). A readout with a
seed_offset is trained with seed k + offset and still written under seed<k>. train_seconds runs from the start of training
to the predicted test set (training, development checks, saving, the reload check, prediction), not reading the inputs.
The paper's runs took 3 to 25 minutes per readout on one H100 (on ArchipelaGO about 18 minutes for the ESM2-15B readout
and 23 for the ProtT5 readout). Host memory: the dense uint8 label matrix and the FP16 embeddings, peak about 26 GB for
archipelago_600k (687,797 proteins) and 9 GB for ArchipelaGO.
"""
import argparse
import json
import time

import numpy as np
import torch

from archipelago.io import read_matrix
from retrain.common import Ontology, add_argument, check_disjoint, data_root, fit_folder, load, log
from retrain.finch import config

NOT_RECIPE = ('encoder', 'embedding', 'seed_offset')  # a readout's entry in configs/finch.json besides its recipe
TOLERANCE = 1e-5


def recipe_of(settings, readout):
    """A readout's training recipe (engine.train's keys): the shared settings and its own."""
    entry = settings['readouts'][readout]
    return {**settings['training'], **{k: v for k, v in entry.items() if k not in NOT_RECIPE}}


def vectors(stem, ids):
    """FP16 embeddings of `ids`, in that order; fails on any missing protein or non-finite vector."""
    _, x = read_matrix(stem, ids)
    if x.dtype != np.float16:
        raise ValueError(f'{stem}: embeddings are {x.dtype}, expected float16')
    out = np.array(x)
    for start in range(0, len(out), 8192):
        if not np.isfinite(out[start:start + 8192]).all():
            raise ValueError(f'{stem}: non-finite vectors')
    return out


class Inputs:
    """One embedding's vectors and one training set's labels on the GPU, with development and test inputs (no test
    labels). allow_shared: training proteins may also be development or test proteins (a published method's own split,
    used as released: retrain.outside); their number is logged."""
    def __init__(self, root, training_set, embedding, allow_shared=False):
        started = time.monotonic()
        self.embedding, self.training_set = embedding, training_set
        train = load(root, training_set, sequences=False)
        dev = load(root, 'development', sequences=False)
        test = load(root, 'test', labels=False, sequences=False)
        self.ontology = Ontology(root)
        self.ontology.assert_closed(train.labels)  # every training positive has its is_a parents
        self.shared = check_disjoint(train, dev, test, allow_shared)
        if train.labels.shape[1] != len(self.ontology.vocabulary):
            raise ValueError(f'{training_set}: {train.labels.shape[1]} label columns, the vocabulary has '
                             f'{len(self.ontology.vocabulary)}')
        from retrain.finch import engine  # CUDA and Triton: imported once the inputs have passed the checks above
        stem = root / 'inputs' / 'embeddings' / embedding
        self.n = len(train.seq_ids)
        self.x = torch.as_tensor(vectors(stem / training_set, train.seq_ids), device='cuda')  # float16, widened per batch
        self.y = torch.empty(train.labels.shape, dtype=torch.uint8, device='cuda')
        for start in range(0, self.n, 32768):
            self.y[start:start + 32768] = torch.as_tensor(train.labels[start:start + 32768].toarray(), device='cuda')
        self.dev_x = torch.as_tensor(vectors(stem / 'development', dev.seq_ids), device='cuda').float()
        self.test_x = torch.as_tensor(vectors(stem / 'test', test.seq_ids), device='cuda').float()
        self.ontology.assert_closed(dev.labels)
        self.dev_y = dev.labels.toarray().astype(np.uint8)
        self.metric = engine.MicroAUPRC(self.dev_y)
        log(f'inputs {embedding} on {training_set}: {self.n:,} training proteins ({self.shared:,} also development or '
            f'test proteins), {len(dev.seq_ids):,} development, {len(test.seq_ids):,} test, '
            f'{torch.cuda.memory_allocated() / 1e9:.1f} GB on the GPU, {time.monotonic() - started:.0f} s')


def fit(root, training_set, seed, readout, recipe, inputs, seed_offset=0):
    """Train one readout with seed + seed_offset, save, reload and reproduce, predict test, write outputs."""
    from retrain.finch import engine
    r = recipe
    trained_seed = seed + seed_offset
    out = fit_folder(root, 'finch', training_set, seed) / readout
    log(f'{readout} on {training_set}, seed {trained_seed}: {inputs.n:,} training proteins, {inputs.x.shape[1]} dimensions, '
        f'{inputs.y.shape[1]} terms')
    started = time.monotonic()

    def check(event, best, best_step):
        log('check', readout, training_set, trained_seed, json.dumps(dict(event, best=best, best_update=best_step)))

    state, dev_scores, summary = engine.train(inputs.x, inputs.y, inputs.dev_x, inputs.dev_y, r, trained_seed, check)
    out.mkdir(parents=True, exist_ok=True)
    torch.save(dict(state_dict=state, readout=readout, seed=trained_seed, recipe=r, input_dim=int(inputs.x.shape[1]),
                    terms=len(inputs.ontology.vocabulary), vocabulary=inputs.ontology.vocabulary,
                    best_step=summary['best_update']), out / 'readout.pt')
    dev = dev_scores.cpu().numpy()
    head = engine.load_head(torch.load(out / 'readout.pt', map_location='cuda', weights_only=True)['state_dict'])
    reproduced = float(np.abs(engine.predict(head, inputs.dev_x).cpu().numpy() - dev).max())
    if not reproduced <= TOLERANCE:  # NaN fails too
        raise ValueError(f'{out}: reloaded checkpoint differs from its development scores by {reproduced}')
    test = engine.predict(head, inputs.test_x).cpu().numpy()
    propagated = inputs.ontology.propagate(dev)
    metrics = dict(dev_raw_micro_auprc=summary['dev_raw_micro_auprc'], dev_propagated_micro_auprc=inputs.metric(propagated),
                   dev_propagated_micro_fmax=inputs.metric.fmax(propagated))
    train_seconds = time.monotonic() - started
    for split, x in (('dev', dev), ('test', test)):
        np.save(out / f'{split}_scores.npy', x.astype(np.float32))
    (out / 'history.json').write_text(json.dumps(dict(
        readout=readout, training_set=training_set, seed=trained_seed, recipe=r, best_step=summary['best_update'],
        steps=summary['updates'], dev_metrics=metrics, reproduced_dev_max_abs=reproduced, train_seconds=train_seconds,
        median_updates_per_second=summary['median_updates_per_second'], checks=summary['history']), indent=1) + '\n',
        encoding='utf-8')
    log(f"wrote {out}: best development check at update {summary['best_update']:,} of {summary['updates']:,}; "
        f"development micro-AUPRC {metrics['dev_propagated_micro_auprc']:.4f} (propagated), {train_seconds / 60:.1f} min")
    del head
    return metrics


def main():
    settings = config()
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--training-set', required=True, help='a dataset name under datasets/')
    parser.add_argument('--seed', type=int, required=True, help=f"the paper's seeds: {settings['seeds']}")
    parser.add_argument('--readouts', nargs='+', choices=list(settings['readouts']), default=list(settings['weights']))
    parser.add_argument('--allow-shared-proteins', action='store_true',
                        help="train although training proteins are also development or test proteins (a published "
                             "method's own split, used as released: retrain.outside)")
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    if not torch.cuda.is_available():
        raise SystemExit('retrain.finch.train needs a CUDA GPU: the training engine is CUDA graphs and Triton kernels')
    loaded = None
    for readout in a.readouts:  # readouts on the same embedding share its inputs, loaded once
        entry = settings['readouts'][readout]
        if loaded is None or loaded.embedding != entry['embedding']:
            loaded = None
            torch.cuda.empty_cache()
            loaded = Inputs(root, a.training_set, entry['embedding'], a.allow_shared_proteins)
        fit(root, a.training_set, a.seed, readout, recipe_of(settings, readout), loaded, entry.get('seed_offset', 0))


if __name__ == '__main__':
    main()
