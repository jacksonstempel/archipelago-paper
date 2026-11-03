"""Rank the candidate pool with the coverage ranking and write the ranking.

    python scripts/rank_coverage.py [--data DIR] [--picks 600000] [--no-filter] [--devices cuda:0,cuda:1,...]

Reads, under the data root (--data or $ARCHIPELAGO_DATA; data/README.md): datasets/hq (proteins and labels),
datasets/candidates (the screened computational candidates, in the pool's row order, and their labels) and the ESM2-15B
layer-48 embeddings inputs/embeddings/esm2_15b_l48/{hq,candidates}. Writes selection/coverage_ranking.tsv, or with
--no-filter (the ranking without the rare-function filter, a control in Fig. 3b) selection/coverage_ranking_unfiltered.tsv:
rank (1 = first pick), seq_id, score (the pick's largest cosine similarity to the training set when picked).
ArchipelaGO takes the first 150,000 picks; the dose sets take 300,000 and 600,000 (the paper ranked 600,000 with the
filter and 500,000 without it).

This is the paper's ranking stage (archipelago.coverage): every candidate is held on the GPUs, split over --devices
(default: every visible GPU) in consecutive parts of rows, about 92 GB of FP16 vectors for the 8,952,628 candidates; the
filter only marks which may be picked. Ties go to the lower candidate row. The paper ran on four H100 GPUs with 512 picks
per pass and checked the first 2,000 picks against one-at-a-time greedy (--reference-picks); its 600,000 picks took 52
seconds after the vectors were loaded (normalizing the candidates' vectors on the CPU takes a few minutes more). Run on the paper's candidate and HQ vectors with four GPUs and batch 512, this
code reproduced both of the paper's rankings exactly (600,000 filtered and 500,000 unfiltered picks, every score). A
different number of GPUs or batch changes the shapes of the FP16 matrix products, and other input vectors change the
similarities, either of which can reorder near-tied candidates: the paper's ranking read HQ vectors from an earlier
extraction than the released training embeddings (it handled X and U differently; every row differs in the last digits).
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
    parser.add_argument("--batch", type=int, default=512, help="picks considered per pass (the paper's: 512)")
    parser.add_argument("--reference-picks", type=int, default=2000,
                        help="first rank this many picks one at a time and require the batched passes to agree")
    parser.add_argument("--devices", help="comma-separated torch devices (default: every visible GPU)")
    args = parser.parse_args()
    root = data_root(args.data)

    hq = load(root, "hq", sequences=False)
    candidates = load(root, "candidates", sequences=False, labels=not args.no_filter)
    if args.no_filter:
        eligible = np.ones(len(candidates.seq_ids), dtype=bool)
    else:
        if candidates.labels.shape[1] != hq.labels.shape[1]:
            raise ValueError("HQ and candidate labels have different vocabularies")
        hq_counts = np.bincount(hq.labels.indices, minlength=hq.labels.shape[1])
        eligible = coverage.rare_function_filter(candidates.labels.indptr, candidates.labels.indices, hq_counts,
                                                 args.threshold)
        print(f"rare-function filter: {eligible.sum()} of {len(candidates.seq_ids)} candidates carry a term held by "
              f"fewer than {args.threshold} HQ proteins", flush=True)

    embeddings = root / "inputs/embeddings/esm2_15b_l48"
    _, hq_vectors = read_matrix(embeddings / "hq", hq.seq_ids)
    _, vectors = read_matrix(embeddings / "candidates", candidates.seq_ids)
    devices = args.devices.split(",") if args.devices else [f"cuda:{i}" for i in range(torch.cuda.device_count())]
    if not devices:
        raise SystemExit("No GPU visible; pass --devices (e.g. cpu, very slow for the full pool)")
    picks, scores, _ = coverage.rank(coverage.unit_vectors(vectors), torch.as_tensor(coverage.unit_vectors(hq_vectors)),
                                     args.picks, eligible=eligible, devices=devices, batch=args.batch,
                                     reference_picks=args.reference_picks)

    out = root / "selection" / ("coverage_ranking_unfiltered.tsv" if args.no_filter else "coverage_ranking.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["rank", "seq_id", "score"])
        writer.writerows((k + 1, candidates.seq_ids[r], f"{s:.6f}") for k, (r, s) in enumerate(zip(picks, scores)))
    print(f"Wrote {out}", flush=True)


if __name__ == "__main__":
    main()
