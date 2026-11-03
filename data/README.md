# The data root

Every command reads and writes under one folder, the data root, given by `--data DIR` or the environment variable
`ARCHIPELAGO_DATA`. Nothing is read from or written to any other place (apart from model weights that the language-model
libraries download into their own caches).

```
<root>/
  release/            the ArchipelaGO dataset as downloaded, and any files obtained on request (below)
  external/           third-party downloads: the published methods' training data, GO-Plus, Swiss-Prot, ...
  datasets/<name>/    every dataset in the working format (below), written by retrain.datasets and the methods' builders
  selection/          coverage rankings computed with scripts/rank_coverage.py
  inputs/             derived inputs: embeddings/<encoder>/<set>, and one folder per method (deepgose, transfew, dpfunc,
                      structures)
  runs/<method>/<training set>/seed<k>/   every fit's predictions
  scoring/            the scorer's bundle (data/) and one configuration per scoring run
  results/<name>/     the scorer's results and the analyses
  figures/            the figures and tables
  cache/              the scorer's prepared inputs
```

## Conventions

- **Protein identifier** (`seq_id`, `protein_id` in the release): the SHA-256 hex digest of the amino-acid sequence in upper
  case. Identical sequences share one identifier; `accessions` lists the UniProt accessions with that sequence.
- **Labels**: GO molecular-function terms, propagated (every `is_a` ancestor within the vocabulary of an annotated term is
  included) and restricted to the 8,739-term vocabulary. One GO release for everything: `go-basic.obo` of the UniProt-GOA
  files of 2026-07-29 (GO release 2026-07-26), for labels, propagation and information accretion (Smin), for every method.
- **Matrices** (embeddings, features): `<stem>.npy` (proteins x columns) with the row identifiers (seq_ids), one per line,
  in `<stem>.ids.txt` (`archipelago.io.read_matrix`). Embeddings are FP16.

## The dataset release (`release/`)

The Hugging Face dataset, as downloaded (Parquet, one folder per table):

| Table | Rows | Columns |
|---|---|---|
| `training/` | the 237,797 ArchipelaGO proteins | `protein_id`, `accessions`, `sequence`, `part` (`hq` or `selected`), `selection_rank` (1 to 150,000 for selected proteins), `go_terms` (propagated labels) |
| `training_annotations/` | their asserted (direct) annotations | `protein_id`, `go_id`, `evidence_class`, `evidence`, `source`, `assigned_by` |
| `benchmark/development`, `benchmark/test` | 2,500 and 5,000 proteins | `protein_id`, `accessions`, `sequence`, `cluster` (homology cluster), `go_terms` |
| `benchmark_experimental/`, `benchmark_masked/`, `benchmark_all_goa/` | the other label versions | as `benchmark`; `benchmark_masked` adds `masked_terms` |
| `benchmark_annotations/` | their asserted annotations, including the accepted review labels | as `training_annotations` |
| `benchmark_reviews/` | every reviewed non-experimental annotation of the benchmark, accepted or not, with its reasoning | `protein_id`, `go_id`, `decision` (`accept` or `do_not_accept`), the reasoning and sources |
| `benchmark_no_hit/` | the 5,000 test proteins | `protein_id`, `no_hit_up_to` (the largest E-value cutoff at which the protein has no training hit) |
| `go_terms/` | the 8,739 vocabulary terms | `go_id`, `name`, `column` (output column), `parents`, `information_accretion`, protein counts |
| `reference_scores/` | the paper's models on every benchmark version | one row per version, model and training set |
| `ontology/go-basic.obo.gz` | the GO release of every label | |
| `fasta/` | the sequences of the training, development and test proteins | |
| `embeddings/<encoder>/{training,development,test}.npy` | the embeddings the paper's Finch was trained on, `esm2_15b_l48` and `prott5_l21` | float16, rows in the order of `training` and `benchmark` |

`python -m retrain.datasets.release` turns these into `datasets/{ontology,hq,archipelago,development,test}`. The scoring with
other test labels (`retrain.scoring.test_labels`) rebuilds each version from the test annotations and reviews, checks it
against the released table, and scores it from a bundle of the release; in the test annotations, the accepted reviews
have the source `literature review (this work)`.

The released embeddings can stand in for `scripts/embed.py` on `archipelago`, `development` and `test`, so that Finch and
Embedding-kNN can be trained on ArchipelaGO without running the language models. This writes them in the
working format (`<root>` is the data root):

```python
import shutil
from pathlib import Path

import pyarrow.parquet as pq

root = Path("<root>")
for encoder in ("esm2_15b_l48", "prott5_l21"):
    (root / "inputs/embeddings" / encoder).mkdir(parents=True, exist_ok=True)
    for name, table, file in (("archipelago", "training/train", "training"), ("development", "benchmark/development",
                              "development"), ("test", "benchmark/test", "test")):
        stem = root / "inputs/embeddings" / encoder / name
        shutil.copyfile(root / "release/embeddings" / encoder / f"{file}.npy", f"{stem}.npy")
        ids = pq.read_table(root / "release" / f"{table}.parquet", columns=["protein_id"])["protein_id"].to_pylist()
        Path(f"{stem}.ids.txt").write_text("\n".join(ids) + "\n", encoding="utf-8")
```

### Files not in the public release

The public release is the dataset above. Rebuilding some of the paper's other training sets exactly also needs files
that are available from the corresponding author on request; the code reads them from these places when present and
otherwise rebuilds an equivalent (it says so in its log):

| File | Used for |
|---|---|
| `release/candidates/{proteins,annotations}.parquet` and `inputs/embeddings/esm2_15b_l48/candidates.{npy,ids.txt}` | the coverage ranking's 8,952,628 candidates and their embeddings (the ranking itself and the Fig. 3a controls) |
| `release/selection/coverage_ranking.tsv`, `coverage_ranking_unfiltered.tsv` (`rank`, `seq_id`, `score`) | the paper's rankings, so that exactly tied candidates come in the paper's order (dose sets, Fig. 3a control) |
| `release/hq_random/{proteins,annotations}.parquet` | the 150,000 random computational proteins of HQ + random |
| `release/hq_support.tsv` (`seq_id`, `go_id`, `evidence`, `assigned_by`, `route`) | the GOA and experimental-label sets exactly (without it they are approximated from the release annotations) |
| `inputs/structures/chosen.tsv` (`seq_id`, `source`, `entry`, `version`) | the paper's AlphaFold Database models for DPFunc |

## Working formats

### `datasets/<name>/`

| File | Content |
|---|---|
| `proteins.tsv` | one row per protein in the dataset's fixed order: `seq_id`, `accessions` (comma-separated), then dataset-specific columns (`cluster` in development and test; `part` and `selection_rank` in ArchipelaGO and the dose sets; `part` and `order_rank` in the control sets) |
| `sequences.fasta` | `>seq_id` and the sequence, same order |
| `annotations.tsv` | asserted annotations before propagation: `seq_id`, `go_id`, `evidence`, `source` |
| `labels.npz` | propagated labels: a CSR uint8 matrix (proteins x 8,739, columns in vocabulary order), `scipy.sparse.save_npz` |

`datasets/ontology/` holds `go-basic.obo` and `vocabulary.json` (the 8,739 terms, a JSON list in column order).

Names: `ontology`, `development`, `test`, `hq`, `archipelago`, `goa`, `hq_explabels`, `hq_random`, `pool_random`,
`filter_random`, `coverage_unfiltered`, `archipelago_300k`, `archipelago_600k`, `candidates`, and the published methods'
original data `original_deepgose`, `original_transfew`, `original_dpfunc`.

### `runs/<method>/<training set>/seed<k>/`

`dev_scores.npy` and `test_scores.npy`: float32, (2,500 or 5,000) x 8,739, rows in development / test order, columns in
vocabulary order; raw scores, which the scorer propagates (DeepGO-SE's are already propagated, as its published inference
does). `record.json` is written last and marks the fit complete. The published methods' records hold `train_seconds`,
the time spent training (DeepGO-SE: in each world's `record.json`); Finch's is in each readout's `history.json`. Further checkpoints of a fit are
`<split>_scores_<name>.npy` (TransFew: `bcebest`, `auprcbest`). Methods: `finch`, `deepgose`, `transfew`, `dpfunc`,
`naive`, `diamond_knn`, `diamond_knn_screen`, `embedding_knn`, and the recipe ablation's models (`ladder_*` and
`finch_prott5_l24`, listed in `configs/finch.json`); the training set of the original data is `original`.
Finch's readouts are in `runs/finch/<training set>/seed<k>/<readout>/`, with `readout.pt` (`state_dict`, `readout`,
`seed`, `recipe`, `input_dim`, `terms`, `vocabulary`, `best_step`), which `scripts/export_package.py` packages.

### `scoring/`

`scoring/data/` is the scorer's bundle: `vocabulary.txt`, `go-basic.obo`, `benchmark/{development,test}.tsv`
(`protein_id`, `accessions`, `cluster`) with `benchmark/{development,test}_labels.tsv` (`protein_id`, `go_id`,
propagated), and `archipelago/{proteins,labels}.tsv` holding the HQ proteins (`part` = `hq`), from which information
accretion is estimated. `scoring/<run>.json` lists the prediction files of one scoring run (format in the docstring of
`archipelago.scorer`).

### `results/<name>/`

Written by `archipelago.scorer` (its `README.md` describes the files: `metrics.tsv`, `curves.tsv`, `replicates.npz`,
`predictions.tsv`, `inputs.json`) plus `contrasts.tsv`, the paired differences of the seed means with their bootstrap
intervals. `figures/__init__.py` lists which results each figure reads.
