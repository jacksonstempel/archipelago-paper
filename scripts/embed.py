"""Embed the proteins of a dataset with one of Finch's encoders.

    python scripts/embed.py --set archipelago --encoder esm2_15b_l48 [--data DIR] [--device cuda]

Encoders: Finch's two (archipelago.embed: esm2_15b_l48, prott5_l21); esm2_3b_l36, the input of the recipe ablation's
first step (Fig. 4a): ESM2-3B's final layer, 36, extracted as Finch extracts ESM2-15B (archipelago.embed.esm2: windows of
at most 1,022 residues overlapping by at least 256, each residue averaged over its windows, mean over the protein); and
prott5_l24, ProtT5's last layer (the encoder output after its closing LayerNorm), pooled as Finch pools layer 21 (masked
mean over the first 1,024 tokens, special tokens included), the input of the ProtT5-last-layer variant.

Reads datasets/<set>/sequences.fasta under the data root (--data or $ARCHIPELAGO_DATA). Writes
inputs/embeddings/<encoder>/<set>.npy (FP16, proteins x dimension, rows in dataset order) and <set>.ids.txt (seq_ids).
Each chunk is written to the .npy file (a memory map) as it is done, so memory does not grow with the dataset; the
.ids.txt file is written last, once every row is.
"""
import argparse
import re
import time
from pathlib import Path

import numpy as np

from archipelago.embed import LOADERS, MODELS, esm2
from retrain.common import load
from retrain.paths import data_root


def prott5_last(batch=16, device="cuda"):
    """ProtT5's last layer, 24: the encoder output after its closing LayerNorm, tokenized and pooled as
    archipelago.embed.prott5 does for its inner layers."""
    import torch
    from transformers import T5EncoderModel, T5Tokenizer
    name = MODELS["prott5_l21"]["model"]
    tokenizer = T5Tokenizer.from_pretrained(name, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(name, torch_dtype=torch.float16).eval().to(device)
    if len(model.encoder.block) != 24:
        raise ValueError(f"Expected 24 encoder blocks, found {len(model.encoder.block)}")

    @torch.inference_mode()
    def embed(sequences):
        out = []
        for a in range(0, len(sequences), batch):
            texts = [" ".join(re.sub(r"[UZOB]", "X", s.upper())) for s in sequences[a:a + batch]]
            encoded = tokenizer(texts, add_special_tokens=True, padding=True, truncation=True, max_length=1024,
                                return_tensors="pt").to(device)
            states = model(**encoded).last_hidden_state
            mask = encoded["attention_mask"].float().unsqueeze(-1)
            out.append(((states.float() * mask).sum(1) / mask.sum(1).clamp_min(1.0)).cpu().numpy())
        return np.concatenate(out)

    return embed


ENCODERS = {**LOADERS,
            "esm2_3b_l36": lambda device: esm2(layer=36, model_name="facebook/esm2_t36_3B_UR50D", device=device),
            "prott5_l24": lambda device: prott5_last(device=device)}


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
