"""How far each Fig. 3b rule's 150,000 added proteins lie from HQ (Fig. 3a).

    python -m retrain.scoring.distances [--data DIR] [--device cuda]

For every added protein of the five training sets of Fig. 3b (hq_random, pool_random, filter_random, coverage_unfiltered,
archipelago: the rows whose `part` is not `hq`), the distance to its nearest HQ protein in ESM2-15B embedding space: one
minus the highest cosine similarity, with every embedding scaled to unit length and the products taken in FP16 on the
GPU, as in the paper. For reference, each HQ protein's distance to its nearest other HQ protein (`hq`).

Reads datasets/<set>/proteins.tsv and inputs/embeddings/esm2_15b_l48/<set>.{npy,ids.txt} (scripts/embed.py) for hq and the
five sets. Writes results/distances/:
  distances.npz   every protein's distance, one array per set
  summary.tsv     per set: proteins, minimum, quartiles, the share within 0.05 and beyond 0.20
  densities.tsv   condition, distance, density: the smoothed density that Fig. 3a draws. The distances are 1 - an FP16
                  cosine, so below 0.5 they come in steps of 2^-11; the histogram bins are 2^-9 wide and aligned at 0 (four
                  steps each; unaligned bins alias into ripples), smoothed with a Gaussian kernel (sigma 0.004) reflected
                  at 0, where distances stop.
Seconds on one GPU.
"""
import argparse

import numpy as np
import torch

from archipelago.io import read_matrix
from retrain.common import add_argument, data_root, load, log

SETS = ('hq_random', 'pool_random', 'filter_random', 'coverage_unfiltered', 'archipelago')
ADDED = 150_000
BIN, SIGMA, TOP = 2.0 ** -9, 0.004, 0.6


def unit(x, device):
    t = torch.as_tensor(np.asarray(x), device=device).float()
    return (t / t.norm(dim=1, keepdim=True)).half()


def nearest(queries, bank, exclude_self=False, chunk=8192):
    """1 - the highest cosine of each row of queries against bank (both unit, FP16); exclude_self: queries is bank, and
    each row's own column is skipped."""
    out = torch.empty(len(queries), device=queries.device)
    for start in range(0, len(queries), chunk):
        similarity = queries[start:start + chunk] @ bank.T
        if exclude_self:
            rows = torch.arange(similarity.shape[0], device=queries.device)
            similarity[rows, rows + start] = -2
        out[start:start + chunk] = 1 - similarity.max(dim=1).values.float()
    return out.clamp_min(0).cpu().numpy()


def density(x):
    """The smoothed density of distances x at the bin centres (module docstring)."""
    edges = np.arange(0, int(TOP / BIN) + 1) * BIN
    half = int(4 * SIGMA / BIN)
    kernel = np.exp(-0.5 * (np.arange(-half, half + 1) * BIN / SIGMA) ** 2)
    kernel /= kernel.sum()
    counts = np.histogram(x, edges)[0].astype(float)
    padded = np.r_[counts[half:0:-1], counts]  # reflected at 0
    return (edges[:-1] + edges[1:]) / 2, np.convolve(padded, kernel, mode='same')[half:] / (len(x) * BIN)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--device', default='cuda')
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    out = root / 'results' / 'distances'
    out.mkdir(parents=True, exist_ok=False)
    embeddings = root / 'inputs' / 'embeddings' / 'esm2_15b_l48'
    hq_ids = load(root, 'hq', labels=False, sequences=False).seq_ids
    hq = unit(read_matrix(embeddings / 'hq', hq_ids)[1], a.device)
    distances = {'hq': nearest(hq, hq, exclude_self=True)}
    for name in SETS:
        ds = load(root, name, labels=False, sequences=False)
        added = [h for h, part in zip(ds.seq_ids, ds.columns['part']) if part != 'hq']
        if len(added) != ADDED or set(added) & set(hq_ids):
            raise ValueError(f'{name}: {len(added)} added proteins (expected {ADDED}), or some of them are HQ proteins')
        distances[name] = nearest(unit(read_matrix(embeddings / name, added)[1], a.device), hq)
        log(f'{name}: median distance to the nearest HQ protein {np.median(distances[name]):.4f}')
    np.savez(out / 'distances.npz', **distances)
    with open(out / 'summary.tsv', 'w', newline='\n') as handle:
        handle.write('condition\tproteins\tmin\tq25\tmedian\tq75\tshare_within_0.05\tshare_beyond_0.20\n')
        for name, d in distances.items():
            q25, median, q75 = np.quantile(d, [0.25, 0.5, 0.75])
            handle.write(f'{name}\t{len(d)}\t{d.min():.4f}\t{q25:.4f}\t{median:.4f}\t{q75:.4f}\t'
                         f'{(d <= 0.05).mean():.4f}\t{(d > 0.20).mean():.4f}\n')
    with open(out / 'densities.tsv', 'w', newline='\n') as handle:
        handle.write('condition\tdistance\tdensity\n')
        for name, d in distances.items():
            centres, values = density(d)
            handle.writelines(f'{name}\t{c:.5f}\t{v:.4f}\n' for c, v in zip(centres, values))
    log(f'wrote {out}')


if __name__ == '__main__':
    main()
