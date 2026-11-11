"""TransFew re-tuned on ArchipelaGO (Supplementary Note 8, Supplementary Table 15): the search's fits, its log, the choice and
the scoring of the chosen setting.

The search keeps the paper's TransFew retraining on ArchipelaGO (retrain.transfew.train: features, label GCN, frequent/rare
split, loss, Adam with weight decay 1e-4, cosine schedule with T_max 50, batch 256, 150 epochs) and varies two settings:
the learning rate (1e-5, 3e-5, 1e-4, 3e-4, 1e-3; published 1e-4) and the dropout after the two hidden layers of the two
protein networks (0.2, 0.4, 0.6, 0.8; published 0.2), one fit per setting with seed 0. Every fit keeps three checkpoints
(retrain.transfew.train); the selection uses the epoch with the highest propagated development micro-AUPRC (`auprcbest`,
the trainer's per-epoch metric, the first epoch reaching the maximum). The chosen setting is the complete seed-0 fit with
the highest such score (ties: the smaller learning rate, then the larger dropout). Had the best value of a setting been at
the edge of its range (dropout 0.2, the published value, excepted), one further step would have been searched; for the
paper's search it was not (the choice, 1e-4 and 0.6, lies inside both ranges). The chosen setting is then trained with
seeds 1 and 2 (seed 0 is the search's fit), and only these three fits are scored on the test set. No test label is read
by the trainer or here before `score`.

    python -m retrain.transfew.retune grid [--seed K]     the commands of the search's fits (or of seed K at the choice)
    python -m retrain.transfew.retune log                 every fit folder, complete or not, and the current best
    python -m retrain.transfew.retune choose              writes retune_transfew/chosen.json once (fails if it exists)
    python -m retrain.transfew.retune score [--device D]  scores the chosen setting's three seeds against TransFew with
                                                          the published recipe (three checkpoint rules) and Finch, all on
                                                          ArchipelaGO, on the main test labels and on the reviewed ones

Fits: <root>/retune_transfew/runs/lr<lr>_d<dropout>/seed<k>/, written by
`python -m retrain.transfew.train archipelago <k> --lr <lr> --dropout <dropout> --out <that folder>` (one GPU, about 38
minutes each on an H100; the paper's 22 fits took 13.9 H100-hours). `log` writes retune_transfew/search_log.tsv.
`score` writes the scorer config scoring/retune.json (the methods transfew_retuned, the chosen fits' best-development-AUPRC
checkpoint, and transfew_retuned_last, the same fits' last epoch, beside finch, transfew, transfew_bcebest and
transfew_auprcbest from runs/), then results_masked/retune (the main test labels: the scorer's --labels main, as
retrain.scoring.score --labels masked) and results/retune (the reviewed labels), each with contrasts.tsv, both from the
released dataset (release/). Needs the main run's fits.
"""
import argparse
import csv
import json
from datetime import datetime, timezone

from retrain.common import add_argument, data_root, log, score_file
from retrain.scoring.score import contrasts, dataset, results_dir, run_scorer, write_config

LRS = ('1e-5', '3e-5', '1e-4', '3e-4', '1e-3')
DROPOUTS = ('0.2', '0.4', '0.6', '0.8')
SEEDS = (0, 1, 2)            # the seeds of every other TransFew fit
COLUMNS = ('lr', 'dropout', 'seed', 'status', 'epochs_done', 'dev_auprc_best', 'best_epoch', 'dev_auprc_last',
           'dev_auprc_lowest_loss', 'lowest_loss_epoch', 'wall_minutes', 'folder')


def folder(root):
    return root / 'retune_transfew'


def fits(root):
    """One row per fit folder. A complete fit has record.json; any other folder is listed from its history.json."""
    rows = []
    for fit in sorted(folder(root).glob('runs/lr*_d*/seed*')):
        lr, dropout = fit.parent.name[2:].split('_d')
        row = dict(lr=float(lr), dropout=float(dropout), seed=int(fit.name[4:]), folder=fit.relative_to(root).as_posix())
        history = json.loads((fit / 'history.json').read_text()) if (fit / 'history.json').exists() else []
        row['epochs_done'] = len(history)
        if (fit / 'record.json').exists():
            r = json.loads((fit / 'record.json').read_text())
            if (r['lr'], r['hidden_dropout'], r['seed'], r['training_set']) != (row['lr'], row['dropout'], row['seed'],
                                                                               'archipelago'):
                raise ValueError(f'{fit}: the record is of another fit')
            best, low = r['side']['auprcbest']['epoch'], r['side']['bcebest']['epoch']
            row.update(status='complete', dev_auprc_best=history[best]['dev_auprc_propagated'], best_epoch=best,
                       dev_auprc_last=history[-1]['dev_auprc_propagated'],
                       dev_auprc_lowest_loss=history[low]['dev_auprc_propagated'], lowest_loss_epoch=low,
                       wall_minutes=r['wall_seconds'] / 60)
        else:
            row['status'] = 'incomplete'
        rows.append(row)
    return rows


def best_of(rows, seed):
    done = [r for r in rows if r['status'] == 'complete' and r['seed'] == seed]
    if not done:
        raise ValueError(f'no complete fit of seed {seed}')
    best = max(done, key=lambda r: (r['dev_auprc_best'], -r['lr'], r['dropout']))
    lrs, dropouts = sorted({r['lr'] for r in done}), sorted({r['dropout'] for r in done})
    edges = dict(lr_lowest=best['lr'] == lrs[0], lr_highest=best['lr'] == lrs[-1],
                 dropout_lowest=best['dropout'] == dropouts[0], dropout_highest=best['dropout'] == dropouts[-1])
    return best, edges, dict(lr=lrs, dropout=dropouts)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('mode', choices=['grid', 'log', 'choose', 'score'])
    parser.add_argument('--seed', type=int, help='grid: the commands of this seed at the chosen setting')
    parser.add_argument('--device', choices=['cuda', 'cpu'], help="score: the scorer's device")
    add_argument(parser)
    a = parser.parse_args(argv)
    root = data_root(a.data)
    if a.mode == 'grid':
        if a.seed is None:
            points, seed = [(lr, d) for lr in LRS for d in DROPOUTS], SEEDS[0]
        else:
            chosen = json.loads((folder(root) / 'chosen.json').read_text())
            points, seed = [(chosen['lr_folder'], chosen['dropout_folder'])], a.seed
        for lr, d in points:
            out = folder(root) / 'runs' / f'lr{lr}_d{d}' / f'seed{seed}'
            print(f'python -m retrain.transfew.train archipelago {seed} --lr {lr} --dropout {d} --out {out} --data {root}')
        return
    rows = fits(root) if a.mode in ('log', 'choose') else None
    if a.mode == 'log':
        path = folder(root) / 'search_log.tsv'
        with open(path, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter='\t', lineterminator='\n', restval='')
            writer.writeheader()
            writer.writerows(rows)
        for r in rows:
            print(*(r.get(c, '') if not isinstance(r.get(c), float) else f'{r[c]:.6g}' for c in COLUMNS[:-1]), sep='\t')
        best, edges, searched = best_of(rows, SEEDS[0])
        print('BEST', json.dumps(dict(lr=best['lr'], dropout=best['dropout'], dev_auprc_best=best['dev_auprc_best'],
                                      edges=edges, searched=searched)))
        log(f'wrote {path}')
    elif a.mode == 'choose':
        if any(r['status'] != 'complete' for r in rows):
            raise ValueError('an incomplete fit folder exists: finish or account for it first')
        best, edges, searched = best_of(rows, SEEDS[0])
        lr_folder, dropout_folder = best['folder'].split('/')[-2][2:].split('_d')
        chosen = dict(lr=best['lr'], dropout=best['dropout'], lr_folder=lr_folder, dropout_folder=dropout_folder,
                      search_seed=SEEDS[0], dev_auprc_best=best['dev_auprc_best'],
                      checkpoint='auprcbest (highest propagated development micro-AUPRC)', edges=edges, searched=searched,
                      fits=[{c: r.get(c) for c in COLUMNS} for r in rows if r['seed'] == SEEDS[0]],
                      written=datetime.now(timezone.utc).isoformat(timespec='seconds'))
        with open(folder(root) / 'chosen.json', 'x') as handle:
            handle.write(json.dumps(chosen, indent=1) + '\n')
        print(json.dumps({k: v for k, v in chosen.items() if k != 'fits'}))
    else:  # score
        chosen = json.loads((folder(root) / 'chosen.json').read_text())
        fit = folder(root) / 'runs' / f"lr{chosen['lr_folder']}_d{chosen['dropout_folder']}"
        models = []
        for method, side in (('finch', None), ('transfew', None), ('transfew', 'bcebest'), ('transfew', 'auprcbest')):
            seeds = sorted((root / 'runs' / method / 'archipelago').glob('seed*'), key=lambda p: int(p.name[4:]))
            models.append((method + (f'_{side}' if side else ''), 'archipelago',
                           {s.name[4:]: {split: score_file(s, prefix, side) for split, prefix in
                                         (('development', 'dev'), ('test', 'test'))} for s in seeds}))
        for method, side in (('transfew_retuned', 'auprcbest'), ('transfew_retuned_last', None)):
            seeds = {}
            for seed in SEEDS:
                if not (fit / f'seed{seed}' / 'record.json').exists():
                    raise FileNotFoundError(f"{fit / f'seed{seed}'}: incomplete fit of the chosen setting")
                seeds[str(seed)] = {split: score_file(fit / f'seed{seed}', prefix, side) for split, prefix in
                                    (('development', 'dev'), ('test', 'test'))}
            models.append((method, 'archipelago', seeds))
        for labels in ('masked', 'reviewed'):
            config = write_config(root / 'scoring' / f'retune_{labels}.json', dataset(root), models)
            out = results_dir(root, labels) / 'retune'
            run_scorer(root, config, out, a.device, labels=labels)
            contrasts(out)


if __name__ == '__main__':
    main()
