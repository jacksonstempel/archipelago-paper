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

