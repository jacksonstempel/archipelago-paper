"""Finch's training engine (retrain.finch.engine): the loss kernel, checkpoint keeping and reloading, and that a run is
reproducible. Needs a CUDA GPU and Triton; skipped without them."""
import numpy as np
import pytest
import torch

try:
    import triton  # noqa: F401
    GPU = torch.cuda.is_available()
except ImportError:
    GPU = False
pytestmark = pytest.mark.skipif(not GPU, reason='the Finch engine needs a CUDA GPU and Triton')

RECIPE = dict(width=32, dropout=0.1, lr=1e-3, rho=0.5, eta=0.01, batch=32, warmup=10, check_interval=20, patience=3,
              weight_decay=0.0, eps=1e-8, loss='bce')


def problem(seed=0):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(200, 16)).astype(np.float16)
    y = (x[:, :4] > 0.5).astype(np.uint8)
    return torch.as_tensor(x[:150], device='cuda'), torch.as_tensor(y[:150], device='cuda'), \
        torch.as_tensor(x[150:], device='cuda').float(), y[150:]


def test_bce_kernel_is_binary_cross_entropy():
    from retrain.finch.engine import loss_gradient
    z = torch.randn(5, 128, device='cuda') * 3  # 3 terms padded to 128 columns: padding never contributes
    y = (torch.rand(5, 3, device='cuda') < 0.4).to(torch.uint8)
    loss, grad = loss_gradient(z, y, 'bce')
    expected = torch.nn.functional.binary_cross_entropy_with_logits(z[:, :3], y.float())
    torch.testing.assert_close(loss, expected, rtol=1e-5, atol=1e-7)
    torch.testing.assert_close(grad[:, :3], (torch.sigmoid(z[:, :3]) - y.float()) / 15, rtol=1e-5, atol=1e-8)
    assert not grad[:, 3:].any()


@pytest.mark.parametrize('width, dropout, rho', [(32, 0.1, 0.5), (None, None, None)], ids=['finch', 'linear, AdamW alone'])
def test_training_keeps_the_best_checkpoint_and_reloads_it(width, dropout, rho):
    from retrain.finch.engine import MicroAUPRC, load_head, predict, train
    x, y, dev_x, dev_y = problem()
    recipe = dict(RECIPE, width=width, dropout=dropout, rho=rho)
    state, scores, summary = train(x, y, dev_x, dev_y, recipe, seed=1)
    best = max(h['dev_raw_micro_auprc'] for h in summary['history'])
    assert summary['best_update'] == next(h['update'] for h in summary['history'] if h['dev_raw_micro_auprc'] == best)
    assert summary['dev_raw_micro_auprc'] == best == MicroAUPRC(dev_y)(scores)
    assert torch.equal(predict(load_head(state), dev_x), scores)


def test_a_run_is_reproducible():
    from retrain.finch.engine import train
    x, y, dev_x, dev_y = problem()
    first, _, a = train(x, y, dev_x, dev_y, RECIPE, seed=7)
    second, _, b = train(x, y, dev_x, dev_y, RECIPE, seed=7)
    assert a['best_update'] == b['best_update'] and all(torch.equal(first[k], second[k]) for k in first)
