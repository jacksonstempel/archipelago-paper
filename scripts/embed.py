"""Embed the proteins of a dataset with one of Finch's encoders.

    python scripts/embed.py --set archipelago --encoder esm2_15b_l48 [--data DIR] [--device cuda]

Encoders: Finch's two (archipelago.embed: esm2_15b_l48, prott5_l21) and esm2_3b_l36, the input of the recipe ablation's
first step (Fig. 4a): ESM2-3B's final layer, 36, extracted as Finch extracts ESM2-15B (archipelago.embed.esm2: windows of
at most 1,022 residues overlapping by at least 256, each residue averaged over its windows, mean over the protein).

Reads datasets/<set>/sequences.fasta under the data root (--data or $ARCHIPELAGO_DATA). Writes
inputs/embeddings/<encoder>/<set>.npy (FP16, proteins x dimension, rows in dataset order) and <set>.ids.txt (seq_ids).
Each chunk is written to the .npy file (a memory map) as it is done, so memory does not grow with the dataset; the
.ids.txt file is written last, once every row is.
"""
import argparse
import time
from pathlib import Path

import numpy as np

from archipelago.embed import LOADERS, esm2
from retrain.common import load
from retrain.paths import data_root

ENCODERS = {**LOADERS,
            "esm2_3b_l36": lambda device: esm2(layer=36, model_name="facebook/esm2_t36_3B_UR50D", device=device)}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--set", required=True, help="a dataset folder name under datasets/")
    parser.add_argument("--encoder", choices=sorted(ENCODERS), required=True)
    parser.add_argument("--data", help="data root (default: $ARCHIPELAGO_DATA)")
    parser.add_argument("--chunk", type=int, default=1000, help="proteins per progress report")
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    root = data_root(args.data)

    dataset = load(root, args.set, labels=False)
    sequences = [s.upper() for s in dataset.sequences]
    out = root / "inputs/embeddings" / args.encoder / args.set
    out.parent.mkdir(parents=True, exist_ok=True)
    Path(f"{out}.ids.txt").unlink(missing_ok=True)  # written last: its presence marks a complete matrix
    embed = ENCODERS[args.encoder](device=args.device)
    stored, started = None, time.monotonic()
    for a in range(0, len(sequences), args.chunk):
        x = np.asarray(embed(sequences[a:a + args.chunk])).astype(np.float16)
        if not np.isfinite(x).all():
            raise ValueError(f"Embeddings of proteins {a} to {a + len(x) - 1} overflow FP16")
        if stored is None:
            stored = np.lib.format.open_memmap(f"{out}.npy", mode="w+", dtype=np.float16,
                                               shape=(len(sequences), x.shape[1]))
        stored[a:a + len(x)] = x
        done = a + len(x)
        rate = done / max(time.monotonic() - started, 1e-6)
        print(f"{args.set} {args.encoder}: {done}/{len(sequences)} proteins, {rate:.1f}/s, "
              f"eta {(len(sequences) - done) / rate / 60:.1f} min", flush=True)
    stored.flush()
    Path(f"{out}.ids.txt").write_text("\n".join(dataset.seq_ids) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
