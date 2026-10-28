"""Train one TransFew fit (training set, seed) with the published recipe and write its predictions.

Recipe (official training.py with hparams.mf['full'], commit 5a4c9d7): 150 epochs; batch 256, shuffled each epoch;
Adam lr 1e-4, weight decay 1e-4 (coupled L2); CosineAnnealingLR(T_max=50) stepped once per epoch (the learning rate
rises again after epoch 50, as in the source); loss = mean of elementwise BCE times get_weights(batch labels, 3); the
label GCN frozen (its dropout and BatchNorm still act in train mode, as in the source). Seeds 0, 1, 2 (the source fixes
17).

Frequent/rare split from this training set's label counts: frequent = more than 30 proteins (the source's
term_indicies[30]; its upper cap of 70,000 only removed the MF root, which the vocabulary lacks), rare = the rest,
including terms without a training positive (the source has no output for those; here they are trained as negatives so
that every vocabulary term gets a score).

Checkpoints. The primary prediction is the model after the last epoch, as the released inference loads it
(predict.py, current_checkpoint.pt; the paper does not say). Two more are kept for the Supplementary Information: the
epoch with the lowest development BCE (sum of the unweighted batch means over shuffled batches of 256, `<=`; the
training script's rule; `bcebest`) and the epoch with the highest propagated development micro-AUPRC (`auprcbest`).
After training, a fresh model loads each checkpoint and must reproduce the development scores saved with it before
its test scores are written.

    python -m retrain.transfew.train TRAINING_SET SEED [--data DIR]

TRAINING_SET: original (datasets/original_transfew), goa, hq or archipelago. Inputs: retrain.transfew.features and
retrain.transfew.label_gcn. Output: runs/transfew/<training set>/seed<k>/ (raw scores, vocabulary order). One GPU;
about 10 min (original data) to 1 h (ArchipelaGO) on an H100.
"""
import argparse
import json
import math
import os
import random
import time

import numpy as np
import torch

from archipelago.io import read_matrix
from retrain import common
from retrain.transfew import features, label_gcn

RECIPE = dict(epochs=150, batch_size=256, lr=1e-4, weight_decay=1e-4, cosine_t_max=50, weight_factor=3,
              frequent_threshold=30)
SIDE = ('bcebest', 'auprcbest')


def seed_all(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False


def micro_ap(scores, labels):
    """Micro-averaged average precision of a dense score matrix (torch), as sklearn's
    average_precision_score(labels.ravel(), scores.ravel()): sum over distinct thresholds of (R_n - R_{n-1}) P_n."""
    s, y = scores.reshape(-1).double(), labels.reshape(-1).double()
    order = torch.argsort(s, descending=True)
    s, y = s[order], y[order]
    tp, fp = torch.cumsum(y, 0), torch.cumsum(1 - y, 0)
    last = torch.ones_like(s, dtype=torch.bool)
    last[:-1] = s[1:] != s[:-1]          # last index of each run of equal scores
    tp, fp = tp[last], fp[last]
    precision, recall = tp / (tp + fp), tp / tp[-1]
    return float(torch.sum(torch.diff(recall, prepend=recall.new_zeros(1)) * precision))


def predict(model, x, batch=256):
    model.eval()
    with torch.no_grad():
        return torch.cat([model(x[i:i + batch]) for i in range(0, len(x), batch)])


def dev_bce(scores, labels, batch=256):
    """The source's validation loss: shuffled batches of `batch`, sum of the batch means of unweighted BCE. In eval mode
    the model scores each protein independently, so batching the scores equals batching the inputs."""
    perm = torch.randperm(len(scores)).to(scores.device)
    return sum(float(torch.nn.functional.binary_cross_entropy(scores[idx], labels[idx]))
               for idx in perm.split(batch))


def split(counts, threshold):
    """Frequent (more than `threshold` training proteins) and rare output columns."""
    freq, rare = np.flatnonzero(counts > threshold), np.flatnonzero(counts <= threshold)
    if not len(freq) or not len(rare):
        raise ValueError('empty frequent or rare group')
    return torch.as_tensor(freq), torch.as_tensor(rare)


def build(root, n_terms, freq, rare):
    from retrain.transfew.model import TFun
    gcn = torch.load(label_gcn.folder(root) / 'gcn.pt', map_location='cpu')
    graph = torch.load(label_gcn.folder(root) / 'graph.pt', map_location='cpu')
    model = TFun(n_terms, freq, rare, gcn['state_dict'], graph)
    for name, p in model.named_parameters():      # training.py: label_embedding frozen
        if name.startswith('label_embedding'):
            p.requires_grad = False
    return model


def keep(model, raw, out, name):
    """Save the model as checkpoint_<name>.pt and the development scores it gave at this epoch."""
    torch.save(model.state_dict(), out / f'checkpoint_{name}.pt.tmp')
    os.replace(out / f'checkpoint_{name}.pt.tmp', out / f'checkpoint_{name}.pt')
    np.save(out / f'saved_dev_{name}.npy', raw)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('training_set', choices=['original', 'goa', 'hq', 'archipelago'])
    parser.add_argument('seed', type=int)
    common.add_argument(parser)
    parser.add_argument('--epochs', type=int, default=RECIPE['epochs'], help='fewer epochs for a quick check only')
    args = parser.parse_args()
    root = common.data_root(args.data)
    from retrain.transfew.model import weighted_bce
    began = time.monotonic()
    out = common.fit_folder(root, 'transfew', args.training_set, args.seed)
    if (out / 'record.json').exists():
        common.log(f'{out} is complete; skipped (remove record.json to retrain)')
        return
    out.mkdir(parents=True, exist_ok=True)
    seed_all(args.seed)
    device = torch.device('cuda')

    ontology = common.Ontology(root)
    n_terms = len(ontology.vocabulary)
    dataset = 'original_transfew' if args.training_set == 'original' else args.training_set
    train = common.load(root, dataset)
    dev = common.load(root, 'development')
    test = common.load(root, 'test', labels=False)
    if set(train.seq_ids) & (set(dev.seq_ids) | set(test.seq_ids)):
        raise ValueError(f'{dataset} shares proteins with the development or test set')
    vectors = features.folder(root) / 'vectors'
    x_train, x_dev, x_test = (read_matrix(vectors, d.seq_ids)[1] for d in (train, dev, test))
    if len(train.seq_ids) % RECIPE['batch_size'] == 1:
        raise ValueError('a final training batch of one protein would break BatchNorm in train mode')
    counts = np.asarray(train.labels.sum(0)).ravel()
    freq, rare = split(counts, RECIPE['frequent_threshold'])
    model = build(root, n_terms, freq, rare).to(device)

    xt = torch.as_tensor(np.asarray(x_train, dtype=np.float32)).to(device)
    yt = torch.from_numpy(train.labels.toarray()).to(device)   # uint8 on the GPU; float per batch
    xd = torch.as_tensor(np.asarray(x_dev, dtype=np.float32)).to(device)
    yd = torch.from_numpy(dev.labels.toarray().astype(np.float32)).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=RECIPE['lr'], weight_decay=RECIPE['weight_decay'])
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, RECIPE['cosine_t_max'])
    common.log(f'transfew {args.training_set} seed {args.seed}: {len(xt):,} proteins, {len(freq):,} frequent and '
               f'{len(rare):,} rare terms')

    best_bce, bce_epoch, best_auprc, auprc_epoch, history = math.inf, None, -math.inf, None, []
    for epoch in range(args.epochs):
        model.train()
        train_loss = 0.0
        for idx in torch.randperm(len(xt), device=device).split(RECIPE['batch_size']):
            optimizer.zero_grad()
            loss = weighted_bce(model(xt[idx]), yt[idx].float(), RECIPE['weight_factor'])
            loss.backward()
            optimizer.step()
            train_loss += loss.item()
        scores = predict(model, xd)
        bce = dev_bce(scores, yd)
        scheduler.step()
        if not (math.isfinite(train_loss) and math.isfinite(bce)):
            raise ValueError(f'non-finite loss at epoch {epoch}')
        raw = scores.cpu().numpy()
        auprc = micro_ap(torch.from_numpy(ontology.propagate(raw)).to(device), yd)
        if bce <= best_bce:
            best_bce, bce_epoch = bce, epoch
            keep(model, raw, out, 'bcebest')
        if auprc > best_auprc:
            best_auprc, auprc_epoch = auprc, epoch
            keep(model, raw, out, 'auprcbest')
        history.append(dict(epoch=epoch, train_loss_sum=train_loss, dev_bce_sum=bce, dev_auprc_propagated=auprc))
        (out / 'history.json').write_text(json.dumps(history) + '\n')
        left = (time.monotonic() - began) / (epoch + 1) * (args.epochs - epoch - 1) / 60
        common.log(f'epoch {epoch}: train loss {train_loss:.2f}, dev BCE {bce:.3f}, dev micro-AUPRC {auprc:.4f}, '
                   f'about {left:.0f} min left')
    keep(model, raw, out, 'last')
    del xt, yt

    xs = torch.as_tensor(np.asarray(x_test, dtype=np.float32)).to(device)
    scores, reproduced = {}, {}
    for name in ('last', *SIDE):
        fresh = build(root, n_terms, freq, rare)
        fresh.load_state_dict(torch.load(out / f'checkpoint_{name}.pt', map_location='cpu'))
        fresh.to(device)
        dev_scores = predict(fresh, xd).cpu().numpy()
        reproduced[name] = common.check_reproduced(np.load(out / f'saved_dev_{name}.npy'), dev_scores)
        scores[name] = dev_scores, predict(fresh, xs).cpu().numpy()
        (out / f'saved_dev_{name}.npy').unlink()
    record = dict(method='transfew', training_set=args.training_set, dataset=dataset, seed=args.seed, recipe=RECIPE,
                  epochs_run=args.epochs, primary='last epoch',
                  side=dict(bcebest=dict(epoch=bce_epoch, rule='lowest development BCE'),
                            auprcbest=dict(epoch=auprc_epoch, rule='highest propagated development micro-AUPRC')),
                  frequent_terms=len(freq), rare_terms=len(rare), reproduced_dev_max_abs=reproduced,
                  wall_seconds=time.monotonic() - began)
    common.write_fit(out, *scores['last'], record, side={name: scores[name] for name in SIDE})
    common.log(f'done: BCE-best epoch {bce_epoch}, AUPRC-best epoch {auprc_epoch}, '
               f'{(time.monotonic() - began) / 60:.0f} min')


if __name__ == '__main__':
    main()
