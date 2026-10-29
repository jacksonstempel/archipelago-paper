"""Rank the candidate pool with the coverage ranking and write the ranking.

    python scripts/rank_coverage.py [--data DIR] [--picks 600000] [--no-filter] [--devices cuda:0,cuda:1,...]

Reads, under the data root (--data or $ARCHIPELAGO_DATA; data/README.md): datasets/hq (proteins and labels),
datasets/candidates (the screened computational candidates and their labels) and the ESM2-15B layer-48 embeddings
inputs/embeddings/esm2_15b_l48/{hq,candidates}. Writes selection/coverage_ranking.tsv, or with --no-filter (the
ranking without the rare-function filter, a control in Fig. 3b) selection/coverage_ranking_unfiltered.tsv: rank
(1 = first pick), seq_id, score (the pick's largest cosine similarity to the training set when picked).
ArchipelaGO takes the first 150,000 picks; the dose sets take 300,000 and 600,000, the paper's run. The candidates are
split over --devices (default: every visible GPU): the filtered pool needs about 54 GB of GPU memory, the unfiltered one
about 92 GB plus working space; the paper's 600,000 picks took 51 seconds on four GPUs.

Ties go to the candidate first in a random order seeded by --order-seed. The paper's run broke ties in another fixed
random order, so exactly tied candidates can come in a different order; the released ranking is the reference.
"""
import argparse
import csv

import numpy as np
import torch

from archipelago import coverage
from archipelago.io import read_matrix
from retrain.common import load
from retrain.paths import data_root


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", help="data root (default: $ARCHIPELAGO_DATA)")
    parser.add_argument("--picks", type=int, default=600000)
    parser.add_argument("--threshold", type=int, default=1000,
                        help="rare-function filter: keep candidates with a term held by fewer HQ proteins")
    parser.add_argument("--no-filter", action="store_true", help="rank every candidate (the Fig. 3b control)")
    parser.add_argument("--batch", type=int, default=512, help="picks considered per pass (same result as one at a time)")
    parser.add_argument("--order-seed", type=int, default=0)
    parser.add_argument("--devices", help="comma-separated torch devices (default: every visible GPU)")
    args = parser.parse_args()
    root = data_root(args.data)

    hq = load(root, "hq", sequences=False)
    candidates = load(root, "candidates", sequences=False, labels=not args.no_filter)
    if args.no_filter:
        kept = np.arange(len(candidates.seq_ids))
    else:
        if candidates.labels.shape[1] != hq.labels.shape[1]:
            raise ValueError("HQ and candidate labels have different vocabularies")
        hq_counts = np.bincount(hq.labels.indices, minlength=hq.labels.shape[1])
        y = candidates.labels
        rows, terms = np.repeat(np.arange(y.shape[0]), np.diff(y.indptr)), y.indices
        kept = np.flatnonzero(coverage.rare_function_filter(rows, terms, len(candidates.seq_ids), hq_counts, args.threshold))
        print(f"rare-function filter: {len(kept)} of {len(candidates.seq_ids)} candidates carry a term held by fewer than "
              f"{args.threshold} HQ proteins", flush=True)

    embeddings = root / "inputs/embeddings/esm2_15b_l48"
    _, hq_vectors = read_matrix(embeddings / "hq", hq.seq_ids)
    _, vectors = read_matrix(embeddings / "candidates", candidates.seq_ids)
    devices = args.devices.split(",") if args.devices else [f"cuda:{i}" for i in range(torch.cuda.device_count())]
    if not devices:
        raise SystemExit("No GPU visible; pass --devices (e.g. cpu, very slow for the full pool)")
    x = coverage.unit_fp16_shards(vectors, kept, devices)
    reference = coverage.unit_fp16(hq_vectors, device=devices[0])
    priority = np.random.default_rng(args.order_seed).permutation(len(candidates.seq_ids))[kept]
    order = np.argsort(np.argsort(priority))  # the same tie order, renumbered 0..len(kept)-1
    picks, scores = coverage.rank(x, reference, args.picks, priority=order, batch=args.batch)

    out = root / "selection" / ("coverage_ranking_unfiltered.tsv" if args.no_filter else "coverage_ranking.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["rank", "seq_id", "score"])
        writer.writerows((k + 1, candidates.seq_ids[kept[r]], f"{s:.6f}") for k, (r, s) in enumerate(zip(picks, scores)))
    print(f"Wrote {out}", flush=True)


if __name__ == "__main__":
    main()
