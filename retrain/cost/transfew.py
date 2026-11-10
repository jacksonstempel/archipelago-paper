"""TransFew: the authors' released training loop (training.py, commit 5a4c9d7, submodel 'full', MF hyperparameters) timed on
ArchipelaGO inputs, with and without its debugging and logging switches.

    python -m retrain.cost.transfew off|released [--data DIR]

Data path as released: TransFewDataset (per-protein tensors in Python lists, one __getitem__ per protein), the default torch
DataLoader (batch 256, shuffle, no workers), default collate on the CPU. Per batch as released: labels to the GPU, three
index_selects, get_weights, zero_grad, model, three output index_selects, weighted BCE, backward, step, loss.item().
  off       the above only
  released  plus CUDA_LAUNCH_BLOCKING=1 (set before CUDA starts), torch.autograd.set_detect_anomaly(True) around the epoch's
            training loop, and compute_scores (sklearn roc_curve + auc on the CPU, four .item() syncs) for the frequent and
            the rare columns after every training and every validation batch.
The model is the paper's port (retrain/transfew/model.py: the released TFun with the ESM branch, same layers, loss and
optimizer), the optimizer Adam lr 1e-4 weight decay 1e-4 with CosineAnnealingLR(50) as released. Inputs: a random subset of
SUBSET ArchipelaGO training proteins (their ESM2-15B vectors and labels, retrain.transfew.features); the frequent/rare split
comes from the full ArchipelaGO label counts. Per-batch time does not depend on the number of proteins, so epoch times are
scaled to 237,797 proteins (retrain.cost.table). Validation: the released validation loop on the DEVN development proteins.
One H100, a few minutes per mode; needs scikit-learn.
"""
import os
import sys
import time

MODE = sys.argv[1] if len(sys.argv) > 1 else None
if MODE == 'released':
    os.environ['CUDA_LAUNCH_BLOCKING'] = '1'  # training.py line 19; must be set before CUDA starts

import argparse  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402

from archipelago.io import read_matrix  # noqa: E402
from retrain.common import add_argument, data_root, load  # noqa: E402
from retrain.transfew import features  # noqa: E402
from retrain.transfew import train as T  # noqa: E402

SUBSET, DEVN, EPOCHS = 15360, 2500, 3


class TransFewDataset(torch.utils.data.Dataset):
    """MyDataset.TransFewDataset for submodel 'full' (ESM only here): per-protein tensors, one item per __getitem__."""
    def __init__(self, x, y):
        self.esm_features = [torch.from_numpy(v) for v in x]
        self.labels = torch.from_numpy(y)

    def __getitem__(self, index):
        return self.esm_features[index], self.labels[index]

    def __len__(self):
        return len(self.labels)


def get_weights(labels, alpha):  # training.py
    label_frequencies = torch.sum(labels, dim=0)
    total_count = torch.sum(label_frequencies, dim=0)
    label_frequencies = torch.add(label_frequencies, 1e-2)
    return torch.pow(torch.log2(torch.div(total_count, label_frequencies)), alpha)


def compute_scores(labels, preds):  # training.py, verbatim
    from sklearn.metrics import auc, roc_curve
    conf = preds / labels
    true_positives = torch.sum(conf == 1).item() + 1e-5
    false_positives = torch.sum(conf == float('inf')).item()
    true_negatives = torch.sum(torch.isnan(conf)).item()
    false_negatives = torch.sum(conf == 0).item()
    accuracy = (true_positives + true_negatives) / (1.0 * (true_positives + true_negatives + false_positives + false_negatives))
    recall = true_positives / (1.0 * (true_positives + false_negatives))
    precision = true_positives / (1.0 * (true_positives + false_positives))
    fscore = 2 * precision * recall / (precision + recall)
    fpr, tpr, _ = roc_curve(labels.flatten().cpu(), preds.flatten().cpu())
    return accuracy, precision, recall, fscore, auc(fpr, tpr)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('mode', choices=['off', 'released'])
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    mode = a.mode
    torch.backends.cuda.matmul.allow_tf32 = False
    device = 'cuda:0'

    def log(*args):
        print(time.strftime('%H:%M:%S'), mode, *args, flush=True)

    train = load(root, 'archipelago', sequences=False)
    dev = load(root, 'development', sequences=False)
    n_terms = train.labels.shape[1]
    counts = np.asarray(train.labels.sum(0)).ravel()
    freq, rare = T.split(counts, 30)
    rng = np.random.default_rng(0)
    rows = np.sort(rng.choice(len(train.seq_ids), size=SUBSET, replace=False))
    vectors = features.folder(root) / 'vectors'
    x = np.ascontiguousarray(read_matrix(vectors, [train.seq_ids[i] for i in rows])[1], dtype=np.float32)
    y = train.labels[rows].toarray().astype(np.float32)
    x_dev = np.ascontiguousarray(read_matrix(vectors, dev.seq_ids[:DEVN])[1], dtype=np.float32)
    y_dev = dev.labels[:DEVN].toarray().astype(np.float32)
    log('subset', SUBSET, 'dev', DEVN, 'freq/rare', len(freq), len(rare), 'gpu', torch.cuda.get_device_name(0))

    torch.manual_seed(0)
    model = T.build(root, n_terms, freq, rare).to(device)
    freq_ix, rare_ix, full_ix = freq.to(device), rare.to(device), torch.arange(n_terms, device=device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-4, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, 50)
    criterion = torch.nn.BCELoss(reduction='none')
    trainloader = torch.utils.data.DataLoader(TransFewDataset(x, y), batch_size=256, shuffle=True)
    valloader = torch.utils.data.DataLoader(TransFewDataset(x_dev, y_dev), batch_size=256, shuffle=True)

    # the data path alone (collate of 256 items), no GPU work
    t = time.time()
    for i, _ in zip(range(20), trainloader):
        pass
    log('DATA PATH alone ms/batch', round(1000 * (time.time() - t) / 20, 1))

    diag = mode == 'released'
    for epoch in range(EPOCHS):
        t0 = time.time()
        model.train()
        with torch.autograd.set_detect_anomaly(diag):
            k = 0
            for feats, labels in trainloader:
                labels = labels.to(device)
                freq_labels = torch.index_select(labels, 1, freq_ix)
                rare_labels = torch.index_select(labels, 1, rare_ix)
                labels = torch.index_select(labels, 1, full_ix)
                class_weights = get_weights(labels, 3)
                optimizer.zero_grad()
                output = model(feats.to(device))
                output_freq = torch.index_select(output, 1, freq_ix)
                output_rare = torch.index_select(output, 1, rare_ix)
                output = torch.index_select(output, 1, full_ix)
                loss = (criterion(output, labels) * class_weights).mean()
                loss.backward()
                optimizer.step()
                loss.data.item()
                if diag:
                    compute_scores(freq_labels, output_freq > 0.5)
                    compute_scores(rare_labels, output_rare > .5)
                k += 1
        torch.cuda.synchronize()
        train_s = time.time() - t0
        t1 = time.time()
        with torch.no_grad():
            model.eval()
            v = 0
            for feats, labels in valloader:
                labels = labels.to(device)
                freq_labels = torch.index_select(labels, 1, freq_ix)
                rare_labels = torch.index_select(labels, 1, rare_ix)
                output = model(feats.to(device))
                output_freq = torch.index_select(output, 1, freq_ix)
                output_rare = torch.index_select(output, 1, rare_ix)
                loss = criterion(output, labels).mean()
                loss.data.item()
                if diag:
                    compute_scores(freq_labels, output_freq > 0.5)
                    compute_scores(rare_labels, output_rare > .5)
                v += 1
        val_s = time.time() - t1
        scheduler.step()
        log(f'EPOCH {epoch} train batches {k} seconds {train_s:.1f} ms/batch {1000 * train_s / k:.1f} | validation batches {v} '
            f'seconds {val_s:.2f} ms/batch {1000 * val_s / v:.1f}')


if __name__ == '__main__':
    main()
