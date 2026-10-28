"""DeepGO-SE (Kulmanov et al. 2024) retrained on one training set: ten "worlds", then their ensemble.

The published recipe (deepgo2 train.py, commit 6c77c98): ESM2-3B input (features.py); DeepGOModel (model.py) with BCE plus
all four EL losses over the GO-Plus normal forms on every batch; Adam, learning rate 5e-4, batch 37, MultiStepLR at epochs
5 and 20 (x0.1), FP32; the checkpoint with the lowest validation loss (the BCE, averaged over batches); ten worlds (seeds
0 to 9) per training set; the six with the lowest validation loss averaged after propagating each one's scores (train.py
propagates its predictions, evaluate_entailment.py averages the top six). Here:
- the development set is the validation set;
- training stops after 20 epochs without a lower development BCE (at most 128 epochs, the official count);
- the training step is engine.py's (the same computation as the official PyTorch step, on Triton kernels);
- the axioms are GO-Plus of the GO release of the labels, normalized by axioms/Normalizer.groovy;
- each world's scores are propagated as the official propagate_annots does (is_a and every `relationship:` line of
  datasets/ontology/go-basic.obo), then averaged;
- every output covers the 8,739 vocabulary terms.

    python -m retrain.deepgose.train world TRAINING_SET WORLD      (one GPU; WORLD 0..9)
    python -m retrain.deepgose.train ensemble TRAINING_SET          (CPU, after the ten worlds)

TRAINING_SET is a dataset name, or `original` for datasets/original_deepgose. Worlds go to
runs/deepgose/<training set>/seed0/members/world<k>/ (best.pt, history.json, raw scores, record.json), the ensemble to
runs/deepgose/<training set>/seed0/ (propagated scores). A world whose record.json exists is complete and is skipped. On
an H100 a world took 0.6 GPU-hours on average (21 to 90 epochs); the paper's 40 worlds took 24 GPU-hours. The world
stage needs Triton, which compiles its kernels against the Python development headers: install them (python3-dev or
equivalent), or point CPATH at the folder holding Python.h.
"""
import argparse
import json
import math
import random
import time
from pathlib import Path

import numpy as np

from retrain.common import (Ontology, add_argument, check_reproduced, data_root, fit_folder, load, log, write_fit)
from retrain.deepgose import features

METHOD = 'deepgose'
WORLDS, ENSEMBLE = 10, 6
MAX_EPOCHS, PATIENCE = 128, 20
BATCH, LR, MILESTONES, GAMMA = 37, 5e-4, [5, 20], 0.1


def dataset_name(training_set):
    return 'original_deepgose' if training_set == 'original' else training_set


def world_folder(root, training_set, k):
    return fit_folder(root, METHOD, training_set, 0) / 'members' / f'world{k}'


def should_stop(epoch, best_epoch, patience=PATIENCE, max_epochs=MAX_EPOCHS):
    """After `epoch` (0-based): stop when `patience` epochs have passed without a lower development BCE, or at the cap."""
    return epoch - best_epoch >= patience or epoch + 1 >= max_epochs


def micro_auprc(scores, labels):
    """Raw micro-AUPRC on the GPU, for the training log only."""
    import torch
    s, y = scores.flatten(), labels.flatten().double()
    y = y[torch.argsort(s, descending=True)]
    precision = torch.cumsum(y, 0) / torch.arange(1, len(y) + 1, device=y.device, dtype=torch.float64)
    return float((precision * y).sum() / y.sum())


def world(root, training_set, k):
    import torch
    from retrain.deepgose.engine import GraphPlan, StreamingStep
    from retrain.deepgose.model import DeepGOModel, load_normal_forms

    folder = world_folder(root, training_set, k)
    if (folder / 'record.json').exists():
        log(f'{folder} is complete')
        return
    folder.mkdir(parents=True, exist_ok=True)
    began = time.monotonic()
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False

    ontology = Ontology(root)
    terms = ontology.index
    train = load(root, dataset_name(training_set))
    dev = load(root, 'development')
    test = load(root, 'test', labels=False)
    if train.labels.shape[1] != len(terms) or dev.labels.shape[1] != len(terms):
        raise ValueError('labels do not have one column per vocabulary term')
    if set(train.seq_ids) & (set(dev.seq_ids) | set(test.seq_ids)):
        raise ValueError(f'{train.name} shares proteins with the development or test set')
    x, dx, tx = (features.read(root, ds) for ds in (train, dev, test))
    y = train.labels.toarray().astype(np.uint8)

    # Seeding, axioms, model, optimizer and data in the order of the official train.py.
    random.seed(k)
    np.random.seed(k)
    torch.manual_seed(k)
    torch.cuda.manual_seed_all(k)
    forms = load_normal_forms(str(root / 'inputs' / 'deepgose' / 'axioms' / 'go-plus.norm'), terms)
    normal = tuple(torch.tensor(nf, dtype=torch.int64, device='cuda') for nf in forms[:4])
    net = DeepGOModel(features.DIM, len(terms), len(forms[5]), len(forms[4]), 'cuda:0').cuda()
    streaming = StreamingStep(net, GraphPlan(net, normal), lr=LR)
    optimizer = streaming.other
    scheduler = torch.optim.lr_scheduler.MultiStepLR(optimizer, milestones=MILESTONES, gamma=GAMMA)
    n = len(x)
    xg, yg = torch.from_numpy(x).cuda(), torch.from_numpy(y).cuda()
    dxg = torch.from_numpy(dx).cuda()
    dyg = torch.from_numpy(dev.labels.toarray().astype(np.float32)).cuda()
    log(f'{training_set} world {k}: {n:,} proteins, {math.ceil(n / BATCH):,} batches per epoch, normal forms '
        f'{[len(f) for f in forms[:4]]}, {len(forms[5]):,} classes outside the vocabulary, {len(forms[4])} relations')

    def predict(model, inputs):
        model.eval()
        with torch.inference_mode():
            return torch.cat([model(inputs[i:i + BATCH]) for i in range(0, len(inputs), BATCH)])

    best, best_epoch, history, captured_lr = float('inf'), -1, [], None
    for epoch in range(MAX_EPOCHS):
        epoch_start = time.monotonic()
        net.train()
        rows = torch.randperm(n).cuda()  # the official loader's shuffle, applied cumulatively to the data order
        xg, yg = xg[rows], yg[rows]
        train_bce = torch.zeros((), device='cuda', dtype=torch.float64)
        train_el = torch.zeros((), device='cuda', dtype=torch.float64)
        for i in range(0, n, BATCH):
            fx, fy = xg[i:i + BATCH], yg[i:i + BATCH].float()
            # The learning rate is a captured kernel argument: recapture after each milestone; a partial batch runs eagerly.
            if captured_lr != optimizer.param_groups[0]['lr']:
                if len(fx) != BATCH:
                    raise ValueError(f'cannot capture the training step: fewer than {BATCH} training proteins')
                streaming.graph, streaming.lr = None, optimizer.param_groups[0]['lr']
                streaming.capture(fx, fy)
                captured_lr = streaming.lr
            loss = streaming(fx, fy) if len(fx) == BATCH else streaming._step(fx, fy)
            torch._assert_async(torch.isfinite(loss), 'non-finite objective')
            train_bce += streaming.last_bce.double()
            train_el += streaming.last_el.double()
        scores = predict(net, dxg)
        dev_bce = float(torch.stack([torch.nn.functional.binary_cross_entropy(scores[i:i + BATCH], dyg[i:i + BATCH]).double()
                                     for i in range(0, len(scores), BATCH)]).mean())
        if not math.isfinite(dev_bce):
            raise ValueError(f'epoch {epoch}: non-finite development scores')
        if dev_bce < best:
            best, best_epoch = dev_bce, epoch
            best_scores = scores.cpu().numpy()
            torch.save({key: value.detach().cpu() for key, value in net.state_dict().items()}, folder / 'best.pt')
        scheduler.step()
        steps = math.ceil(n / BATCH)
        history.append(dict(epoch=epoch, train_bce=float(train_bce) / steps, train_el=float(train_el) / steps,
                            dev_bce=dev_bce, dev_raw_micro_auprc=micro_auprc(scores, dyg), best_epoch=best_epoch,
                            seconds=time.monotonic() - epoch_start))
        (folder / 'history.json').write_text(json.dumps(history, indent=1) + '\n')
        log(f'epoch {epoch}: dev BCE {dev_bce:.5f} (best {best:.5f} at {best_epoch}), dev raw micro-AUPRC '
            f"{history[-1]['dev_raw_micro_auprc']:.4f}, {history[-1]['seconds']:.0f} s; at most "
            f"{max(0, best_epoch + PATIENCE - epoch) * history[-1]['seconds'] / 60:.0f} min until the stop if no improvement")
        if should_stop(epoch, best_epoch):
            break

    # The kept checkpoint, reloaded, must reproduce the development scores saved when it was selected.
    del streaming, optimizer, scheduler, xg, yg
    torch.cuda.empty_cache()
    fresh = DeepGOModel(features.DIM, len(terms), len(forms[5]), len(forms[4]), 'cuda:0').cuda()
    fresh.load_state_dict(torch.load(folder / 'best.pt', map_location='cuda'))
    reproduced = check_reproduced(best_scores, predict(fresh, dxg).cpu().numpy())
    test_scores = predict(fresh, torch.from_numpy(tx).cuda()).cpu().numpy()
    write_fit(folder, best_scores, test_scores, dict(
        method=METHOD, training_set=training_set, world=k, propagated=False, best_dev_bce=best, selected_epoch=best_epoch,
        epochs_run=epoch + 1, training_proteins=n, reproduced_dev_max_abs=reproduced, gpu=torch.cuda.get_device_name(),
        wall_seconds=time.monotonic() - began,
        selection=f'lowest development BCE: epoch {best_epoch} (0-based) of {epoch + 1}; stop after {PATIENCE} epochs '
                  f'without improvement, at most {MAX_EPOCHS}'))
    log(f'world {k} of {training_set} done: epoch {best_epoch} kept of {epoch + 1}, '
        f'{(time.monotonic() - began) / 60:.0f} min')


def official_ancestors(obo):
    """deepgo2 deepgo/utils.py Ontology(filename, with_rels=True): load() and get_ancestors(), as there: every `is_a` and
    every `relationship:` line (any type) is a parent; alt_ids point to their term; obsolete terms are removed. Returns
    get_ancestors (the term itself included)."""
    ont, obj = dict(), None
    with open(obo, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line == '[Term]':
                if obj is not None:
                    ont[obj['id']] = obj
                obj = dict(is_a=list(), alt_ids=list(), is_obsolete=False)
                continue
            elif line == '[Typedef]':
                if obj is not None:
                    ont[obj['id']] = obj
                obj = None
            else:
                if obj is None:
                    continue
                item = line.split(': ')
                if item[0] == 'id':
                    obj['id'] = item[1]
                elif item[0] == 'alt_id':
                    obj['alt_ids'].append(item[1])
                elif item[0] == 'is_a':
                    obj['is_a'].append(item[1].split(' ! ')[0])
                elif item[0] == 'relationship':
                    obj['is_a'].append(item[1].split()[1])  # add all types of relationships
                elif item[0] == 'is_obsolete' and item[1] == 'true':
                    obj['is_obsolete'] = True
        if obj is not None:
            ont[obj['id']] = obj
    for term_id in list(ont.keys()):
        for t_id in ont[term_id]['alt_ids']:
            ont[t_id] = ont[term_id]
        if ont[term_id]['is_obsolete']:
            del ont[term_id]
    cache = {}

    def get_ancestors(term_id):
        if term_id not in ont:
            return set()
        if term_id not in cache:
            term_set, q = set(), [term_id]
            while q:
                t_id = q.pop()
                if t_id not in term_set:
                    term_set.add(t_id)
                    q += [p for p in ont[t_id]['is_a'] if p in ont]
            cache[term_id] = term_set
        return cache[term_id]
    return get_ancestors


class OfficialPropagation:
    """The official propagate_annots (deepgo/utils.py) for a whole score matrix: each vocabulary term gets the maximum score
    over the vocabulary terms that have it among their ancestors (is_a and every relationship, through the full ontology,
    itself included)."""

    def __init__(self, obo, vocabulary):
        get_ancestors = official_ancestors(obo)
        index = {t: i for i, t in enumerate(vocabulary)}
        below = [[] for _ in vocabulary]
        for t, j in index.items():
            ancestors = get_ancestors(t)
            if t not in ancestors:
                raise ValueError(f'{t} is missing from {obo}')
            for a in ancestors:
                if a in index:
                    below[index[a]].append(j)
        self.below = [np.array(sorted(b)) for b in below]

    def __call__(self, scores):
        scores = np.asfortranarray(scores, dtype=np.float32)
        out = np.empty_like(scores)
        for a, b in enumerate(self.below):
            out[:, a] = scores[:, b].max(axis=1)
        return np.ascontiguousarray(out)


def select_worlds(records, size=ENSEMBLE):
    """The `size` worlds with the lowest development BCE (official evaluate_entailment.get_top_models); ties by index."""
    return [r['world'] for r in sorted(records, key=lambda r: (r['best_dev_bce'], r['world']))[:size]]


def ensemble(root, training_set):
    records = []
    for k in range(WORLDS):
        path = world_folder(root, training_set, k) / 'record.json'
        if not path.exists():
            raise FileNotFoundError(f'world {k} of {training_set} is not complete ({path})')
        records.append(json.loads(path.read_text()))
    chosen = select_worlds(records)
    ontology = Ontology(root)
    propagate = OfficialPropagation(ontology.obo, ontology.vocabulary)
    out = {}
    for split in ('dev', 'test'):
        total = sum(propagate(np.load(world_folder(root, training_set, k) / f'{split}_scores.npy')).astype(np.float64)
                    for k in chosen)
        out[split] = (total / len(chosen)).astype(np.float32)
    write_fit(fit_folder(root, METHOD, training_set, 0), out['dev'], out['test'], dict(
        method=METHOD, training_set=training_set, seed=0, propagated=True, ensemble_worlds=chosen,
        worlds={r['world']: dict(best_dev_bce=r['best_dev_bce'], selected_epoch=r['selected_epoch'],
                                 epochs_run=r['epochs_run']) for r in records},
        selection=f'the {ENSEMBLE} of {WORLDS} worlds with the lowest development BCE, each propagated (is_a and every '
                  'relationship of go-basic.obo, as deepgo2 propagate_annots), then averaged'))
    log(f'ensemble of {training_set}: worlds {chosen}')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('stage', choices=['world', 'ensemble'])
    parser.add_argument('training_set', help='dataset name, or `original` for original_deepgose')
    parser.add_argument('world', nargs='?', type=int, choices=range(WORLDS), help='world (seed) 0-9; stage world only')
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    if a.stage == 'world':
        if a.world is None:
            parser.error('world needs a world number')
        world(root, a.training_set, a.world)
    else:
        ensemble(root, a.training_set)


if __name__ == '__main__':
    main()
