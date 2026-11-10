"""TransFew retraining: the pieces whose errors would silently change results (CPU, tiny inputs)."""
import numpy as np
import pytest
import torch

from retrain.transfew import train
from retrain.transfew.model import get_weights


def test_frequent_and_rare_columns_cover_every_term_once():
    counts = np.array([0, 31, 30, 100, 5])
    freq, rare = train.split(counts, 30)
    assert freq.tolist() == [1, 3] and rare.tolist() == [0, 2, 4]  # frequent means strictly more than 30


def test_loss_weights_follow_the_source_formula():
    """Mirrors TransFew training.py get_weights (commit 5a4c9d7): log2(total batch label count / (column count + 0.01)) ** 3."""
    labels = torch.tensor([[1., 0., 1.], [1., 0., 0.]])
    expected = torch.log2(3 / (torch.tensor([2., 0., 1.]) + 0.01)) ** 3
    assert torch.allclose(get_weights(labels, 3), expected)


def test_micro_ap_equals_sklearn_with_ties():
    metrics = pytest.importorskip('sklearn.metrics')
    rng = np.random.default_rng(0)
    scores = rng.integers(0, 5, size=(40, 6)) / 4  # many ties
    labels = (rng.random((40, 6)) < 0.3).astype(np.float64)
    expected = metrics.average_precision_score(labels.ravel(), scores.ravel())
    assert train.micro_ap(torch.tensor(scores), torch.tensor(labels)) == pytest.approx(expected, abs=1e-12)


def test_the_hidden_dropout_acts_after_both_hidden_layers_only():
    from retrain.transfew.model import TFun_submodel
    published = [m.p for m in TFun_submodel(7).modules() if isinstance(m, torch.nn.Dropout)]
    retuned = [m.p for m in TFun_submodel(7, hidden_dropout=0.6).modules() if isinstance(m, torch.nn.Dropout)]
    assert published == [0.1, 0.2, 0.2] and retuned == [0.1, 0.6, 0.6]  # the input dropout stays 0.1


def fake_fit(folder, lr, dropout, seed, dev_auprc, best, low, complete=True):
    """A fit folder as retrain.transfew.train writes it: history.json, and record.json once complete."""
    import json
    folder.mkdir(parents=True)
    history = [dict(epoch=e, dev_auprc_propagated=dev_auprc * (0.9 if e != best else 1.0)) for e in range(5)]
    (folder / 'history.json').write_text(json.dumps(history))
    if complete:
        (folder / 'record.json').write_text(json.dumps(dict(
            lr=lr, hidden_dropout=dropout, seed=seed, training_set='archipelago', wall_seconds=60.0,
            side=dict(auprcbest=dict(epoch=best), bcebest=dict(epoch=low)))))


def test_retune_chooses_the_best_development_auprc_with_the_tie_rules(tmp_path):
    import json
    from retrain.transfew import retune
    runs = tmp_path / 'retune_transfew' / 'runs'
    fake_fit(runs / 'lr1e-4_d0.2' / 'seed0', 1e-4, 0.2, 0, 0.53, 3, 1)
    fake_fit(runs / 'lr1e-4_d0.6' / 'seed0', 1e-4, 0.6, 0, 0.54, 4, 2)
    fake_fit(runs / 'lr3e-5_d0.6' / 'seed0', 3e-5, 0.6, 0, 0.54, 2, 2)  # a tie: the smaller learning rate wins
    fake_fit(runs / 'lr1e-3_d0.8' / 'seed0', 1e-3, 0.8, 0, 0.10, 0, 0)
    rows = retune.fits(tmp_path)
    assert {(r['lr'], r['dropout']): r['dev_auprc_best'] for r in rows}[(1e-4, 0.6)] == pytest.approx(0.54)
    assert {(r['lr'], r['dropout']): r['dev_auprc_last'] for r in rows}[(1e-4, 0.6)] == pytest.approx(0.54)  # best = last
    retune.main(['choose', '--data', str(tmp_path)])
    chosen = json.loads((tmp_path / 'retune_transfew' / 'chosen.json').read_text())
    assert (chosen['lr'], chosen['dropout'], chosen['lr_folder'], chosen['dropout_folder']) == (3e-5, 0.6, '3e-5', '0.6')
    with pytest.raises(FileExistsError):  # the choice is written once
        retune.main(['choose', '--data', str(tmp_path)])


def test_retune_refuses_to_choose_with_an_incomplete_fit(tmp_path):
    from retrain.transfew import retune
    runs = tmp_path / 'retune_transfew' / 'runs'
    fake_fit(runs / 'lr1e-4_d0.2' / 'seed0', 1e-4, 0.2, 0, 0.53, 3, 1)
    fake_fit(runs / 'lr1e-4_d0.4' / 'seed0', 1e-4, 0.4, 0, 0.54, 4, 2, complete=False)
    with pytest.raises(ValueError, match='incomplete'):
        retune.main(['choose', '--data', str(tmp_path)])
