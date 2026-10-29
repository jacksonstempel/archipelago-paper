"""Finch's training recipe (retrain.finch.recipe): the loss, checkpoint keeping, and one update written out by hand."""
import numpy as np
import pytest
import torch

from archipelago import readout
from retrain.finch import recipe as finch_recipe


def test_loss_is_binary_cross_entropy():
    z = torch.tensor([[2.0, -1.0, 0.3]])
    y = torch.tensor([[1, 0, 0]])
    p = torch.sigmoid(z)
    expected = torch.stack([-torch.log(p[0, 0]), -torch.log(1 - p[0, 1]), -torch.log(1 - p[0, 2])]).mean()
    assert torch.allclose(finch_recipe.bce(z, y), expected)


def test_training_runs_and_keeps_the_best_checkpoint():
    rng = np.random.default_rng(0)
    x = rng.normal(size=(200, 16)).astype(np.float16)
    y = (x[:, :4] > 0.5).astype(np.uint8)
    recipe = dict(width=32, dropout=0.1, lr=1e-3, asam_rho=0.5, batch=32, warmup=10, check_interval=20, patience=3)
    result = finch_recipe.train(x[:150], y[:150], x[150:], y[150:], recipe, seed=1, max_updates=100, log=lambda _: None)
    best = max(h["dev_raw_micro_auprc"] for h in result["history"])
    assert result["best_step"] == next(h["step"] for h in result["history"] if h["dev_raw_micro_auprc"] == best)
    model = readout.Readout(16, 32, 4, 0.1)
    model.load_state_dict(result["state"])
    np.testing.assert_array_equal(readout.predict(model, x[150:]), result["dev_scores"])


@pytest.mark.parametrize("width, rho", [(16, 2.0), (None, None)], ids=["finch", "linear readout without ASAM"])
def test_one_update_is_the_recipe_written_out(width, rho):
    """Binary cross-entropy; ASAM (gradient at w + rho s^2 g / ||s g||, s = |w| + eta, applied to w) or, without it, the
    gradient at w; gradient clipping at norm 1; AdamW's first step (-lr g / (|g| + eps)) and the linear warm-up, computed
    by hand on the batch the sampler draws."""
    rng = np.random.default_rng(0)
    x = rng.normal(size=(64, 10)).astype(np.float32)
    y = (rng.random((64, 5)) < 0.3).astype(np.uint8)
    recipe = dict(width=width, dropout=None if width is None else 0.0, lr=1e-2, asam_rho=rho, batch=8, warmup=4,
                  check_interval=1, patience=5)
    seed = 3
    result = finch_recipe.train(x, y, x[:20], y[:20], recipe, seed=seed, max_updates=1, log=lambda _: None)

    torch.manual_seed(seed)
    model = finch_recipe.readout(10, 5, width, recipe["dropout"])  # the same initial weights as train()
    rows = np.searchsorted(np.cumsum(np.full(64, 1 / 64)), np.random.default_rng(seed).random(8), side="right")
    bx, by = torch.as_tensor(x[rows]), torch.as_tensor(y[rows]).float()
    w = [p.detach().clone() for p in model.parameters()]

    def gradient(weights):
        for p, v in zip(model.parameters(), weights):
            p.data.copy_(v)
        model.zero_grad()
        q = torch.sigmoid(model(bx))
        (-(by * torch.log(q) + (1 - by) * torch.log(1 - q))).mean().backward()
        return [p.grad.clone() for p in model.parameters()]

    g = gradient(w)
    if rho is not None:
        s = [v.abs() + 0.01 for v in w]
        norm = torch.sqrt(sum(((gi * si) ** 2).sum() for gi, si in zip(g, s)))
        g = gradient([v + rho * si * si * gi / norm for v, si, gi in zip(w, s, g)])
    total = torch.sqrt(sum((gi ** 2).sum() for gi in g))
    g = [gi * min(1.0, 1.0 / (float(total) + 1e-6)) for gi in g]
    expected = [v - 1e-2 / 4 * gi / (gi.abs() + 1e-8) for v, gi in zip(w, g)]
    for e, name in zip(expected, dict(model.named_parameters())):
        torch.testing.assert_close(result["state"][name], e, rtol=1e-5, atol=1e-6)
