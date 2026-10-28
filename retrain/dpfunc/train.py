"""Train DPFunc on one training set with one seed, then predict the development and test sets.

    python -m retrain.dpfunc.train --training-set {original,goa,hq,archipelago,...} --seed {0,1,2}

The published recipe, as the released DPFunc_main.py runs it (CSUBioGroup/DPFunc at c5faf4df): combine_inter_model with
1,280 hidden units and 4 heads, focal loss (gamma 2), AdamW (learning rate 1e-4, PyTorch's default weight decay 0.01),
batches of 64 shuffled proteins with the last short batch kept, 15 epochs, no schedule, nothing selected on the
development set; the prediction is the mean of the sigmoid scores of the last three epochs (13 to 15). FP32, TF32 off.
Then the published post-processing: every term gets the maximum score of its descendants over is_a and part_of (the GO
release of the data root, within the vocabulary).

Adaptations: every protein is scored on the 8,739-term vocabulary; GraphConv accepts residues without a contact (the
released code would fail on them); seeds 0, 1 and 2 (the released code sets none).

Writes runs/dpfunc/<training set>/seed<k>/: the three epochs' checkpoints and development scores, then (after the
reloaded checkpoints reproduce the saved development scores) dev_scores.npy, test_scores.npy (post-processed),
{dev,test}_scores_raw.npy (the three-epoch means before post-processing) and record.json. Running the
command again after training finished only predicts. One GPU; 23 H100-hours for the paper's 12 fits (the ArchipelaGO
fits are the longest), with up to about 130 GB of host memory for the memory-mapped ArchipelaGO bank. Needs DGL.
"""
import argparse
import json
import mmap
import random
import time

import numpy as np
import scipy.sparse as sp
import torch

from archipelago.evaluate import micro_auprc
from retrain.common import Ontology, add_argument, check_reproduced, data_root, fit_folder, load, log, write_fit
from retrain.dpfunc.common import DIMENSION, N_INTERPRO, SEEDS, TRAINING_SETS, bank_of, dataset_of, inputs

EPOCHS, AVERAGED, BATCH = 15, 3, 64


class Bank:
    """Residue features and contact edges of a bank, memory-mapped."""

    def __init__(self, root, name):
        folder = inputs(root) / 'banks' / name
        self.seq_ids = json.loads((folder / 'seq_ids.json').read_text())
        self.offsets = np.load(folder / 'feature_offsets.npy')
        self.edge_offsets = np.load(folder / 'edge_offsets.npy')
        self.features = np.load(folder / 'features.npy', mmap_mode='r')
        self.edges = np.load(folder / 'edges.npy', mmap_mode='r')
        if self.features.shape != (self.offsets[-1], DIMENSION) or self.edges.shape != (2, self.edge_offsets[-1]):
            raise ValueError(f'{folder}: the arrays do not match their offsets')

    def release(self):
        """Drop mapped pages after a batch, so host memory stays bounded (Linux)."""
        if hasattr(mmap, 'MADV_DONTNEED'):
            self.features._mmap.madvise(mmap.MADV_DONTNEED)
            self.edges._mmap.madvise(mmap.MADV_DONTNEED)


class Graphs(torch.utils.data.Dataset):
    def __init__(self, bank, rows, labels):
        self.bank, self.rows, self.labels = bank, rows, labels

    def __len__(self):
        return len(self.rows)

    def __getitem__(self, i):
        import dgl
        r = int(self.rows[i])
        a, b = int(self.bank.offsets[r]), int(self.bank.offsets[r + 1])
        e0, e1 = int(self.bank.edge_offsets[r]), int(self.bank.edge_offsets[r + 1])
        uv = torch.from_numpy(np.array(self.bank.edges[:, e0:e1]))
        graph = dgl.graph((uv[0], uv[1]), num_nodes=b - a)
        graph.ndata['x'] = torch.from_numpy(np.array(self.bank.features[a:b]))
        return graph, i, self.labels[i].toarray().ravel().astype(np.float32)


class Split:
    """A dataset with its bank rows, domain counts and labels (zero labels when they are not needed)."""

    def __init__(self, root, name, terms, labels=True):
        self.dataset = load(root, name, labels=labels)
        folder = inputs(root) / 'steps' / name
        self.rows = np.load(folder / 'rows.npy')
        self.domains = sp.load_npz(folder / 'domains.npz').tocsr()
        self.bank = Bank(root, bank_of(name))
        n = len(self.dataset.seq_ids)
        if len(self.rows) != n or self.domains.shape != (n, N_INTERPRO):
            raise ValueError(f'{name}: the step inputs do not have one row per protein (run inputs step {name} again)')
        if [self.bank.seq_ids[r] for r in self.rows] != self.dataset.seq_ids:
            raise ValueError(f"{name}: the bank rows do not hold the dataset's proteins in its order")
        self.labels = self.dataset.labels if labels else sp.csr_matrix((n, terms), dtype=np.uint8)

    def loader(self, shuffle=False):
        from dgl.dataloading import GraphDataLoader
        data = Graphs(self.bank, self.rows, self.labels)
        return GraphDataLoader(data, batch_size=BATCH, drop_last=False, shuffle=shuffle)


def seed_everything(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False


def new_model(terms):
    from retrain.dpfunc.model import combine_inter_model
    model = combine_inter_model(inter_size=N_INTERPRO, inter_hid=1280, graph_size=1280, graph_hid=1280, label_num=terms,
                                head=4).cuda()
    model.GNN.conv1.set_allow_zero_in_degree(True)
    model.GNN.conv2.set_allow_zero_in_degree(True)
    return model


def inter_batch(domains, rows):
    """A batch's domain counts as EmbeddingBag input (indices, offsets, weights)."""
    selected = domains[rows]
    return (torch.from_numpy(selected.indices).long().cuda(), torch.from_numpy(selected.indptr).long().cuda(),
            torch.from_numpy(selected.data).float().cuda())


def evaluate(model, split):
    """Sigmoid scores (in dataset order) and the mean focal loss."""
    from retrain.dpfunc.model import FocalLoss
    model.eval()
    loss_fn, scores, total, n = FocalLoss(), [], 0., 0
    with torch.no_grad():
        for graph, rows, y in split.loader():
            graph = graph.to('cuda')
            split.bank.release()
            logits = model(inter_batch(split.domains, rows.numpy()), graph, graph.ndata['x'])
            total += loss_fn(logits, y.cuda()).item() * len(y)
            n += len(y)
            scores.append(torch.sigmoid(logits).cpu().numpy())
    out = np.concatenate(scores)
    if not np.isfinite(out).all():
        raise FloatingPointError('non-finite scores')
    return out, total / n


def published_postprocessing(ontology, scores):
    """DPFunc's post-processing (paper, Methods; DPFunc/evaluation.py: Ontology(with_rels=True) adds part_of to is_a, and
    every ancestor gets the maximum of its descendants' scores), restricted to the vocabulary's columns."""
    parents, current = {}, None
    for line in open(ontology.obo, encoding='utf-8-sig'):
        line = line.rstrip('\n')
        if line.startswith('['):
            current = None
        elif line.startswith('id: '):
            current = line[4:].strip()
            parents[current] = []
        elif current and line.startswith('is_a: '):
            parents[current].append(line[6:].split(' ! ')[0].strip())
        elif current and line.startswith('relationship: part_of '):
            parents[current].append(line.split()[2])
    ancestors = {}

    def up(term):
        if term not in ancestors:
            ancestors[term] = set()
            for p in parents.get(term, []):
                ancestors[term] |= {p} | up(p)
        return ancestors[term]

    raw = np.asarray(scores, dtype=np.float32)
    out = raw.copy()
    for j, term in enumerate(ontology.vocabulary):
        for a in up(term):
            k = ontology.index.get(a)
            if k is not None:
                np.maximum(out[:, k], raw[:, j], out=out[:, k])
    if not np.array_equal(ontology.propagate(out), out):
        raise RuntimeError('post-processed scores are not closed under is_a')
    return out


def train(root, ontology, training_set, seed, folder):
    from retrain.dpfunc.model import FocalLoss
    terms = len(ontology.vocabulary)
    train_split = Split(root, dataset_of(training_set), terms)
    dev = Split(root, 'development', terms)
    ontology.assert_closed(train_split.labels)
    seed_everything(seed)
    loader = train_split.loader(shuffle=True)
    model = new_model(terms)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-4)
    loss_fn, history, kept, started = FocalLoss(), [], [], time.monotonic()
    log(f'DPFunc on {training_set}, seed {seed}: {len(loader.dataset):,} proteins, {EPOCHS} epochs')
    for epoch in range(1, EPOCHS + 1):
        model.train()
        total, n = 0., 0
        for graph, rows, y in loader:
            graph, y = graph.to('cuda'), y.cuda()
            train_split.bank.release()
            loss = loss_fn(model(inter_batch(train_split.domains, rows.numpy()), graph, graph.ndata['x']), y)
            if not torch.isfinite(loss):
                raise FloatingPointError(f'non-finite focal loss in epoch {epoch}')
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            total += loss.item() * len(y)
            n += len(y)
        scores, dev_loss = evaluate(model, dev)
        history.append(dict(epoch=epoch, train_focal=total / n, dev_focal=dev_loss,
                            dev_raw_micro_auprc=micro_auprc(scores, dev.labels.toarray())))
        elapsed = time.monotonic() - started
        log(json.dumps(history[-1]), f'ETA {(EPOCHS - epoch) * elapsed / epoch / 3600:.2f} h')
        if epoch > EPOCHS - AVERAGED:
            torch.save({'epoch': epoch, 'model_state_dict': model.state_dict()}, folder / f'epoch{epoch}.pt')
            kept.append(scores)
    np.save(folder / 'dev_scores_training.npy', np.mean(np.stack(kept), axis=0).astype(np.float32))
    (folder / 'history.json').write_text(json.dumps(dict(history=history, train_seconds=time.monotonic() - started)))
    (folder / 'TRAINED').touch()


def predict(root, ontology, training_set, seed, folder):
    terms = len(ontology.vocabulary)
    dev, test = Split(root, 'development', terms), Split(root, 'test', terms, labels=False)
    dev_parts, test_parts = [], []
    for epoch in range(EPOCHS - AVERAGED + 1, EPOCHS + 1):
        model = new_model(terms)
        model.load_state_dict(torch.load(folder / f'epoch{epoch}.pt', map_location='cuda')['model_state_dict'])
        dev_parts.append(evaluate(model, dev)[0])
        test_parts.append(evaluate(model, test)[0])
    dev_scores = np.mean(np.stack(dev_parts), axis=0).astype(np.float32)
    difference = check_reproduced(np.load(folder / 'dev_scores_training.npy'), dev_scores)
    test_scores = np.mean(np.stack(test_parts), axis=0).astype(np.float32)
    history = json.loads((folder / 'history.json').read_text())
    record = dict(method='dpfunc', training_set=training_set, seed=seed, epochs=EPOCHS,
                  selection=f'mean of the sigmoid scores of epochs {EPOCHS - AVERAGED + 1}-{EPOCHS} (the published rule)',
                  propagated=True, postprocessing='max over is_a and part_of descendants (published)',
                  reproduced_dev_max_abs=difference, history=history['history'], train_seconds=history['train_seconds'])
    write_fit(folder, published_postprocessing(ontology, dev_scores), published_postprocessing(ontology, test_scores),
              record, side={'raw': (dev_scores, test_scores)})  # the 3-epoch means before post-processing
    log(f'wrote {folder}')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--training-set', required=True, choices=TRAINING_SETS)
    parser.add_argument('--seed', type=int, required=True, choices=SEEDS)
    add_argument(parser)
    arguments = parser.parse_args()
    root = data_root(arguments.data)
    torch.set_num_threads(1)
    ontology = Ontology(root)
    folder = fit_folder(root, 'dpfunc', arguments.training_set, arguments.seed)
    folder.mkdir(parents=True, exist_ok=True)
    if not (folder / 'TRAINED').exists():
        train(root, ontology, arguments.training_set, arguments.seed, folder)
    predict(root, ontology, arguments.training_set, arguments.seed, folder)


if __name__ == '__main__':
    main()
