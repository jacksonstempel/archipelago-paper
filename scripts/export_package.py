"""Build Finch's model package, the folder that Finch.from_pretrained loads (format: archipelago/pretrained.py).

    python scripts/export_package.py --out finch_package [--data DIR]

Inputs:
- --config (default configs/finch.json): the readouts, their weights and the runs to package (release_seeds; the
  released Finch is one run, seed 9301).
- --checkpoints: where each readout's checkpoint is, a pattern with {root}, {readout} and {seed}; default
  "{root}/runs/finch/archipelago/seed{seed}/{readout}/readout.pt", as retrain.finch writes them (state_dict, readout,
  seed, recipe, input_dim, terms, vocabulary): one per readout and packaged seed. Each checkpoint's vocabulary (its
  output columns' GO terms, in order) must equal the dataset vocabulary.
- the data root (--data or $ARCHIPELAGO_DATA): datasets/ontology/{vocabulary.json, go-basic.obo} (is_a edges and term
  names), datasets/development (proteins and labels) and each readout's development embeddings,
  inputs/embeddings/<readout>/development.{npy,ids.txt}.

The threshold: the packaged model (archipelago/pretrained.py) scores the development proteins from their embeddings, and
the threshold of its development micro-Fmax (0.01 to 1.00, the first maximum) becomes the package's default cut, as the
paper's F1 chooses its threshold. The test set is not read.

Upload to the Hub afterwards, e.g. `hf upload jstempel/finch finch_package .` (huggingface_hub's command line).
"""
import argparse
import json
from pathlib import Path

import numpy as np
import torch
from safetensors.torch import save_file

from archipelago import __version__, embed
from archipelago.evaluate import micro_fmax
from archipelago.go import read_terms
from archipelago.io import read_matrix, read_vocabulary
from retrain.common import load
from retrain.paths import data_root
from archipelago.pretrained import CONFIG, NAMES, PARENTS, VOCABULARY, Finch, readout_file

ROOT = Path(__file__).resolve().parents[1]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", help="data root (default: $ARCHIPELAGO_DATA)")
    parser.add_argument("--checkpoints", default="{root}/runs/finch/archipelago/seed{seed}/{readout}/readout.pt")
    parser.add_argument("--config", type=Path, default=ROOT / "configs/finch.json")
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    root = data_root(args.data)
    ontology = root / "datasets/ontology"

    config = json.loads(args.config.read_text(encoding="utf-8"))
    weights, seeds = config["weights"], config["release_seeds"]
    vocabulary = read_vocabulary(ontology / "vocabulary.json")
    obo = ontology / "go-basic.obo"
    terms = read_terms(obo)
    missing = [t for t in vocabulary if t not in terms]
    if missing:
        raise ValueError(f"{obo}: {len(missing)} vocabulary terms missing, e.g. {missing[0]}")
    if args.out.exists() and any(args.out.iterdir()):
        raise FileExistsError(f"{args.out} is not empty")
    (args.out / "readouts").mkdir(parents=True, exist_ok=True)

    readouts = {}
    for name in weights:
        recipe = config["readouts"][name]
        dims = set()
        for seed in seeds:
            path = Path(args.checkpoints.format(root=root, readout=name, seed=seed))
            saved = torch.load(path, map_location="cpu", weights_only=True)
            found = (saved["readout"], saved["seed"], saved["recipe"]["width"], saved["recipe"]["dropout"])
            if found != (name, seed, recipe["width"], recipe["dropout"]):
                raise ValueError(f"{path}: readout, seed, width, dropout {found} do not match the config")
            if saved.get("vocabulary") != vocabulary:
                raise ValueError(f"{path}: its output columns are not the GO terms of {ontology / 'vocabulary.json'} "
                                 "(checkpoint vocabulary missing or different)")
            dims.add(saved["input_dim"])
            save_file({k: v.contiguous() for k, v in saved["state_dict"].items()}, str(args.out / readout_file(name, seed)))
            print(f"{name} seed {seed}: {path}", flush=True)
        if len(dims) != 1:
            raise ValueError(f"{name}: seeds differ in input dimension {dims}")
        source = embed.MODELS[name]
        readouts[name] = dict(encoder=recipe["encoder"], language_model=source["model"], layer=source["layer"],
                              input_dim=dims.pop(), width=recipe["width"], dropout=recipe["dropout"])

    column = {t: i for i, t in enumerate(vocabulary)}
    (args.out / VOCABULARY).write_text("\n".join(vocabulary) + "\n", encoding="utf-8")
    edges = [f"{t}\t{vocabulary[j]}" for t in vocabulary for j in sorted({column[p] for p in terms[t][1] if p in column})]
    (args.out / PARENTS).write_text("\n".join(["go_id\tparent_id"] + edges) + "\n", encoding="utf-8")
    (args.out / NAMES).write_text("\n".join(["go_id\tname"] + [f"{t}\t{terms[t][0]}" for t in vocabulary]) + "\n",
                                  encoding="utf-8")
    header = obo.read_text(encoding="utf-8").split("[Term]", 1)[0].splitlines()
    go_release = next((line[len("data-version:"):].strip() for line in header if line.startswith("data-version:")), None)
    package = dict(
        model="Finch",
        description="Finch: the weighted mean of two readouts' raw scores (one-hidden-layer MLPs on mean-pooled "
                    "protein-language-model embeddings), propagated up the GO is_a graph (with several seeds, the mean "
                    "of the seeds' propagated scores). Molecular-function terms only.",
        weights=weights, seeds=seeds, readouts=readouts, terms=len(vocabulary), threshold=None,
        versions=dict(archipelago=__version__, torch=torch.__version__, go_release=go_release))

    def write(values):
        (args.out / CONFIG).write_text(json.dumps(values, indent=1) + "\n", encoding="utf-8")

    write(package)
    development = load(root, "development", sequences=False)
    dev_ids, dev_y = development.seq_ids, development.labels.toarray()
    embeddings = {name: np.asarray(read_matrix(root / "inputs/embeddings" / name / "development", dev_ids)[1],
                                   dtype=np.float32) for name in weights}
    fmax, threshold = micro_fmax(Finch(args.out, args.device).scores_from_embeddings(embeddings), dev_y)
    package.update(threshold=threshold, threshold_selection=f"development micro-Fmax ({fmax:.4f}) of these scores, "
                   f"thresholds 0.01 to 1.00, first maximum; {len(dev_ids)} development proteins")
    write(package)
    print(f"Wrote {args.out}: development micro-Fmax {fmax:.4f} at threshold {threshold:.2f}", flush=True)


if __name__ == "__main__":
    main()
