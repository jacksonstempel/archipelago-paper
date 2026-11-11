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
  scoring/            one scorer configuration per scoring run
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

The Hugging Face dataset (v2.0), as downloaded (Parquet, one folder per table). Every table ends with a `license` column
(`cc-by-4.0`, or `cc-by-nc-4.0` for rows drawn from OrganiZymeDB, KEGG or P450Rdb); the pipeline reads columns by name and
ignores it.

| Table | Rows | Columns |
|---|---|---|
| `proteins/{train,development,test,screened_out}` | the 237,797 ArchipelaGO proteins, 2,500 development and 5,000 test proteins, and the 4,063 HQ proteins the homology screen removed (do not train on them if you evaluate on this benchmark) | `protein_id`, `accessions`, `sequence`, `part` (`hq` or `selected`) and `selection_rank` (train, screened_out), `cluster` (homology cluster) and `no_hit_up_to` (test: the largest E-value cutoff at which the protein has no training hit), `related_*` (screened_out), `go_terms` (propagated labels) with `evidence` (each term's strongest support: `experimental`, `curated`, `computational`, `review_direct`, `review_inferred`, `review_rejected`) |
| `annotations/{train,development,test,screened_out}` | the asserted annotations, and every reviewed GOA assertion of a benchmark protein with the review's decision | `protein_id`, `go_id`, `level`, `source`, `evidence`, `assigned_by`, `decision` (`accept` or `do_not_accept`), the reasoning and sources |
| `go_terms/` | the 8,739 vocabulary terms | `go_id`, `name`, `column` (output column), `parents`, `information_accretion`, protein counts |
| `reference_scores/` | the paper's models on every version of the test labels | one row per version, model and training set |
| `ontology/go-basic.obo.gz` | the GO release of every label | |
| `embeddings/<encoder>/{train,development,test,screened_out}.npy` | the embeddings the paper's Finch was trained on, `esm2_15b_l48` and `prott5_l21` | float16, rows in the order of `proteins` |
| `coverage_ranking/` | the first 600,000 ranks of the paper's coverage ranking (of the 5,287,663 candidates that pass the rare-function filter) | `rank`, `protein_id`, `accessions`, `coverage_score` |
| `control_sets/{hq_random,pool_random,filter_random,coverage_unfiltered}` | the 150,000 proteins each control set of Fig. 3a adds to HQ, in order | `position`, `protein_id`, `accessions`, `order_rank` |

`python -m retrain.datasets.release` turns these into `datasets/{ontology,hq,archipelago,development,test}`; the labels of
the development and test datasets are the reviewed ones (every level but `review_rejected`). The scorer reads `release/`
directly. The scoring with other test labels (`retrain.scoring.test_labels`) rebuilds each version from the test
annotations and the review's decisions, checks it against the released labels of the same version, and scores it; in the
test annotations, the accepted reviews have the source `literature review (this work)`.

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
    for name, split in (("archipelago", "train"), ("development", "development"), ("test", "test")):
        stem = root / "inputs/embeddings" / encoder / name
        shutil.copyfile(root / "release/embeddings" / encoder / f"{split}.npy", f"{stem}.npy")
        ids = pq.read_table(root / "release/proteins" / f"{split}.parquet", columns=["protein_id"])["protein_id"].to_pylist()
        Path(f"{stem}.ids.txt").write_text("\n".join(ids) + "\n", encoding="utf-8")
```

### Deposited ranking and control sets, and files available on request

The release deposits the paper's coverage ranking (its first 600,000 ranks: `coverage_ranking/`) and the
protein identifiers of the four control sets of Fig. 3a (`control_sets/`). `retrain.datasets.selected` reads the deposited
ranking in place of `release/selection/coverage_ranking.tsv`, and every control set and HQ + random that the pipeline
builds is checked against its deposited table (the same proteins in the same order, with the same ranks).

Building those sets, and a few others exactly, also needs files that are available from the corresponding author on
request (the screened candidate pool is large, and the structures are DPFunc inputs); the code reads them from these
places when present and otherwise rebuilds an equivalent (it says so in its log):

| File | Used for |
|---|---|
| `release/candidates/{proteins,annotations}.parquet` and `inputs/embeddings/esm2_15b_l48/candidates.{npy,ids.txt}` | the screened candidate pool: the coverage ranking's 8,952,628 candidates and their embeddings (recomputing the ranking, and the sequences and labels of the dose sets and the Fig. 3a controls) |
| `release/selection/coverage_ranking_unfiltered.tsv` (`rank`, `seq_id`, `score`) | the unfiltered ranking beyond the 150,000 deposited proteins of `coverage_unfiltered` (the control set screens the first 165,000); a ranking recomputed from other HQ or candidate vectors than the ones the paper's ranking read can order near-tied candidates differently |
| `release/hq_random/{proteins,annotations}.parquet` | the sequences and annotations of the 150,000 random computational proteins of HQ + random (their identifiers and accessions are deposited in `control_sets/hq_random`) |
| `release/hq_support.tsv` (`seq_id`, `go_id`, `evidence`, `assigned_by`, `route`) | the GOA and experimental-label sets exactly (without it they are approximated from the release annotations) |
| `inputs/structures/chosen.tsv` (`seq_id`, `source`, `entry`, `version`) and the structures | the paper's AlphaFold Database models and ColabFold structures for DPFunc |

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

`scoring/<run>.json` lists the prediction files of one scoring run and the dataset (`release/`) the scorer reads (format in
the docstring of `archipelago.scorer`).

### `results/<name>/`

Written by `archipelago.scorer` (its `README.md` describes the files: `metrics.tsv`, `curves.tsv`, `replicates.npz`,
`predictions.tsv`, `inputs.json`) plus `contrasts.tsv`, the paired differences of the seed means with their bootstrap
intervals. `figures/__init__.py` lists which results each figure reads.
