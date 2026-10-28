"""Training of Finch's readouts (archipelago.readout.Readout), as in the paper:

- batches of 256 training proteins drawn uniformly with replacement (NumPy generator seeded with the fit's seed);
- asymmetric loss (gamma+ 1, gamma- 4, probability shift 0.05), averaged over proteins and terms;
- AdamW (betas 0.9 / 0.999, eps 1e-8, weight decay 0), learning rate warmed up linearly over the first 1,000 updates and
  then constant, gradients clipped to global norm 1;
- adaptive sharpness-aware minimization (ASAM, elementwise, eta 0.01): the gradient at the weights perturbed by
  rho * s^2 * g / ||s * g|| with s = |w| + eta, with the same dropout mask in both passes, updates the unperturbed weights;
- every 1,000 updates the development set is scored and the raw (unpropagated) micro-AUPRC checked; training stops after
  20 checks without improvement and keeps the best checkpoint.

The paper's runs used fused GPU kernels with BF16 linear layers and FP32 master weights; here a CUDA run uses BF16 autocast
for the training forward pass, and prediction is FP32.
"""
import time

import numpy as np
import torch
import torch.nn.functional as F

from archipelago.evaluate import micro_auprc
from archipelago.readout import Readout, predict

DEFAULTS = dict(batch=256, warmup=1000, check_interval=1000, patience=20, weight_decay=0.0, adam_eps=1e-8,
                max_grad_norm=1.0, gamma_pos=1.0, gamma_neg=4.0, clip=0.05, asam_eta=0.01)


def asymmetric_loss(logits, targets, gamma_pos=1.0, gamma_neg=4.0, clip=0.05):
    """Asymmetric loss (Ridnik et al., 2021), mean over all protein-term pairs."""
    logits = logits.float()
    targets = targets.bool()
    p = torch.sigmoid(logits)
    q = (torch.sigmoid(-logits) + clip).clamp(max=1.0)  # shifted probability of "negative"
    positive = -((1 - p) ** gamma_pos) * F.logsigmoid(logits)
    negative = -((1 - q) ** gamma_neg) * torch.log(q.clamp(min=1e-8))
    return torch.where(targets, positive, negative).mean()


def train(train_x, train_y, dev_x, dev_y, recipe, seed, device="cpu", max_updates=None, log=print):
    """Train one readout. Returns dict(state, dev_scores, best_step, steps, history).

    train_x: embeddings (proteins x dim, any float dtype); train_y: binary labels (proteins x terms);
    dev_x, dev_y: the development set used for early stopping; recipe: width, dropout, lr, asam_rho and optional
    overrides of DEFAULTS. max_updates stops early (for quick checks only).
    """
    r = dict(DEFAULTS, **recipe)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    rng = np.random.default_rng(seed)
    n, dim = train_x.shape
    terms = train_y.shape[1]
    model = Readout(dim, r["width"], terms, r["dropout"]).to(device)
    parameters = list(model.parameters())
    optimizer = torch.optim.AdamW(parameters, lr=r["lr"], betas=(0.9, 0.999), eps=r["adam_eps"],
                                  weight_decay=r["weight_decay"])
    x = torch.from_numpy(np.array(train_x)).to(device)  # a copy: the input may be a read-only memory map
    y = torch.from_numpy(np.array(train_y, dtype=np.uint8)).to(device)
    dev_y = np.asarray(dev_y)
    cdf = np.cumsum(np.full(n, 1.0 / n))
    cdf /= cdf[-1]  # uniform sampling, drawn as in the paper's sampler
    autocast = torch.autocast("cuda", dtype=torch.bfloat16) if torch.device(device).type == "cuda" else \
        torch.autocast("cpu", enabled=False)

    def loss_at(bx, by, mask):
        with autocast:
            logits = model(bx, mask)
        return asymmetric_loss(logits, by, r["gamma_pos"], r["gamma_neg"], r["clip"])

    best, best_step, best_state, best_scores, stale, step, history = -1.0, 0, None, None, 0, 0, []
    started = time.monotonic()
    while stale < r["patience"] and (max_updates is None or step < max_updates):
        if step % r["check_interval"] == 0:
            drawn = np.searchsorted(cdf, rng.random(r["check_interval"] * r["batch"]), side="right")
            drawn = drawn.reshape(r["check_interval"], r["batch"])
        rows = torch.as_tensor(drawn[step % r["check_interval"]], device=device)
        step += 1
        bx, by = x[rows].float(), y[rows]
        mask = F.dropout(torch.ones(len(rows), r["width"], device=device), p=r["dropout"], training=True)

        # ASAM: gradient at the adaptively perturbed weights, same dropout mask, then restore and update
        optimizer.zero_grad(set_to_none=True)
        loss_at(bx, by, mask).backward()
        with torch.no_grad():
            scales = [p.abs() + r["asam_eta"] for p in parameters]
            norm = torch.sqrt(sum(((p.grad * s) ** 2).sum() for p, s in zip(parameters, scales)))
            steps = [r["asam_rho"] * s * s * p.grad / (norm + 1e-12) for p, s in zip(parameters, scales)]
            for p, e in zip(parameters, steps):
                p.add_(e)
        optimizer.zero_grad(set_to_none=True)
        loss = loss_at(bx, by, mask)
        loss.backward()
        with torch.no_grad():
            for p, e in zip(parameters, steps):
                p.sub_(e)
        torch.nn.utils.clip_grad_norm_(parameters, r["max_grad_norm"])
        for group in optimizer.param_groups:
            group["lr"] = r["lr"] * min(step / r["warmup"], 1.0)
        optimizer.step()

        if step % r["check_interval"] and not (max_updates is not None and step == max_updates):
            continue
        scores = predict(model, dev_x, device=device)
        auprc = micro_auprc(scores, dev_y)
        history.append(dict(step=step, dev_raw_micro_auprc=auprc, loss=float(loss.detach())))
        if auprc > best:
            best, best_step, stale, best_scores = auprc, step, 0, scores
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        else:
            stale += 1
        elapsed = time.monotonic() - started
        log(f"step {step}: dev raw micro-AUPRC {auprc:.4f} (best {best:.4f} at {best_step}, {stale} stale checks), "
            f"{step / elapsed:.0f} updates/s")
    return dict(state=best_state, dev_scores=best_scores, best_step=best_step, steps=step, history=history)
