# ArchipelaGO: reproducing the paper

The code behind every number, figure and table of **ArchipelaGO: better training data for protein function prediction**
(Jackson C. Stempel, Jason Andrew Amburn and Jeremy C. Smith; manuscript in preparation): building the training sets,
retraining Finch, DeepGO-SE, TransFew, DPFunc and three baselines on each of them, scoring, and drawing the figures.

**To use the dataset, the Finch predictor or the benchmark scorer, you do not need this repository:** see
[archipelago](https://github.com/jacksonstempel/archipelago), which this repository builds on.

Contact: Jackson C. Stempel, jstempel@vols.utk.edu.

## Install

Python 3.10 or later. Install PyTorch for your hardware first (a CUDA build for GPUs), then, in a checkout:

```bash
pip install -e .                 # the pipeline (installs the archipelago package)
pip install -e .[embed,figures]  # embeddings and the figure scripts; add deepgose, transfew or dpfunc per published method
pip install -e .[test]           # the tests: python -m pytest
```

Install the repository in place (`-e`): the pipeline reads `configs/finch.json` from the checkout. The extras `deepgose`,
`transfew` and `dpfunc` pin the versions the paper's runs used where they were recorded. External programs (DIAMOND,
InterProScan, ColabFold, Java with Groovy) are called from the command line and are found on `PATH`.

## The data root

Everything the pipeline reads and writes lives under one folder, given by `--data DIR` on every command or by the
environment variable `ARCHIPELAGO_DATA`. [data/README.md](data/README.md) describes its layout and every file format.
The starting point is the ArchipelaGO dataset as downloaded from Hugging Face, in `release/`.

## Reproduce

The stages below rebuild every number of the paper from the released dataset and the published methods' releases. Each
command takes `--data DIR` (or reads `ARCHIPELAGO_DATA`) and prints its options with `--help`. Every stage runs as a
plain command; `examples/slurm/stage.sbatch` is an example of running one as a Slurm job. Times are for the paper's
hardware (NVIDIA H100 80 GB GPUs; CPU stages on 16 to 64 cores).

The paper's computations were run with this code's predecessors on a cluster; the pipeline here is the same computation,
cleaned for release and checked on small inputs, not rerun end to end. Finch's readouts were trained with fused GPU
kernels (BF16 linear layers, FP32 master weights) implementing the recipe that `retrain/finch/recipe.py` implements in
plain PyTorch, so retrained numbers match the paper up to run-to-run noise, not bit for bit. The same holds for the
published methods (GPU nondeterminism) and for the structure bank (below).

### 1. Datasets

```bash
python -m retrain.datasets.release          # release/ (the Hugging Face dataset) -> datasets/{ontology,hq,archipelago,development,test}
python -m retrain.datasets.hq_subsets goa   # GOA: HQ proteins with GOA experimental or curator support, their GOA labels
python -m retrain.datasets.hq_subsets hq_explabels   # HQ with only its experimental labels (Supplementary Information)
```

The release conversion stops unless every set's annotations, propagated on the release's GO file, give exactly the
released labels. CPU, minutes.

### 2. The coverage ranking and the sets built from computational proteins

These stages need the candidate pool (8,952,628 proteins; `python -m retrain.datasets.candidates` converts its release
files into `datasets/candidates`) and its ESM2-15B embeddings, which are not part of the public release (data/README.md). The candidates are one random member of each 50%-identity cluster of the
77.9 million TrEMBL proteins with computational molecular-function annotations only (DIAMOND 2.2.8 `cluster --id 50
--mutual-cover 80 --evalue 0.001`), with relatives of the benchmark and the 102 proteins with labels outside the
vocabulary removed.

```bash
python scripts/rank_coverage.py --picks 600000                 # selection/coverage_ranking.tsv (51 s on four GPUs)
python scripts/rank_coverage.py --picks 500000 --no-filter     # selection/coverage_ranking_unfiltered.tsv (Fig. 4b control)
python -m retrain.datasets.selected check-archipelago          # the ranking's first 150,000 = the released selected proteins
python -m retrain.datasets.selected dose                       # archipelago_300k, archipelago_600k
python -m retrain.datasets.selected controls                   # pool_random, filter_random, coverage_unfiltered (DIAMOND)
python -m retrain.datasets.random_draw                         # hq_random, from its released protein list
```

Exactly tied candidates may come in another order than in the paper's run; the dataset builders use the paper's
rankings instead when they are present (data/README.md).

### 3. Finch and the baselines

```bash
python scripts/embed.py --set <set> --encoder esm2_15b_l48     # every training set, development and test (GPU)
python scripts/embed.py --set <set> --encoder prott5_l21
python -m retrain.finch.train --training-set archipelago --seed 9301    # each training set, seeds 9301-9303 (GPU)
python -m retrain.finch.combine --training-set archipelago               # the 0.60 / 0.40 weighted mean
python -m retrain.baselines.naive --training-sets goa hq archipelago
python -m retrain.baselines.diamond_knn --training-sets goa hq archipelago        # DIAMOND
python -m retrain.baselines.embedding_knn --training-sets goa hq archipelago
python scripts/export_package.py --out finch_package           # the model package Finch.from_pretrained loads
```

Embedding: ESM2-15B about 77 GPU-hours for the 9.9 million candidates, a few GPU-hours for the training sets; ProtT5 is
faster. Finch: 5 to 30 minutes per readout on one H100 (about 2 GPU-hours for every training set and seed of the main
comparison). Baselines: minutes.

### 4. The published methods

Each method keeps its published recipe and is trained on its original data (rebuilt from the authors' release and
screened against the benchmark), GOA, HQ and ArchipelaGO. Every fit writes raw development and test scores; a reloaded
checkpoint must reproduce its saved development scores before it predicts the test set.

**DeepGO-SE** (extra `deepgose`; Triton needs the Python development headers; Java 11 and Groovy 2.4 for the axioms)

```bash
python -m retrain.deepgose.original                 # datasets/original_deepgose (DIAMOND)
bash retrain/deepgose/axioms/build.sh               # GO-Plus axioms, normalized (about 80 GB of Java heap)
python -m retrain.deepgose.features development test goa hq archipelago original_deepgose   # ESM2-3B (GPU)
python -m retrain.deepgose.train world <training set> <k>      # k = 0..9, one GPU each (about 0.6 GPU-hours)
python -m retrain.deepgose.train ensemble <training set>       # the six best worlds, propagated and averaged
```

**TransFew** (extra `transfew`)

```bash
python -m retrain.transfew.original                 # datasets/original_transfew (DIAMOND; about 40 GB of RAM)
python -m retrain.transfew.features embed --shard K --shards N && python -m retrain.transfew.features merge   # ESM2-15B (GPU)
python -m retrain.transfew.label_gcn                # label graph and its GCN (GPU, minutes)
python -m retrain.transfew.train <training set> <seed>         # seeds 0-2; 10 minutes to 1 hour per fit
```

**DPFunc** (extra `dpfunc`; InterProScan 5.78-109.0; ColabFold 1.6.3 with MMseqs2 release 18 and its databases for the structures)

```bash
python -m retrain.dpfunc.original                   # datasets/original_dpfunc (DIAMOND)
python -m retrain.structures.afdb requests && python -m retrain.structures.afdb search && python -m retrain.structures.afdb fetch
python -m retrain.structures.colabfold msa --databases DIR && python -m retrain.structures.colabfold fold
python -m retrain.structures.assemble
python -m retrain.dpfunc.embed                      # ESM-1b layer 31 per residue (GPU, about 4 GPU-hours)
python -m retrain.dpfunc.domains scan --interproscan DIR && python -m retrain.dpfunc.domains merge   # about 100 CPU-hours
python -m retrain.dpfunc.inputs graphs              # then: inputs bank <original|archipelago|evaluation>, inputs step <set>
python -m retrain.dpfunc.train --training-set <training set> --seed <seed>   # seeds 0-2; 23 GPU-hours for all 12 fits
```

Structures: an AlphaFold Database model of the exact sequence where one exists (84% of the paper's proteins),
otherwise an AlphaFold2 structure predicted with ColabFold (alphafold2_ptm, model 1, seed 0, 3 recycles, no templates).
The paper's bank reused earlier ColabFold runs, most of them folded with faster kernels (97.3% 3Di identity to standard
folds); a bank built from scratch therefore differs slightly. `afdb fetch --chosen` downloads exactly the paper's
AlphaFold Database entries when the paper's choice list is present (data/README.md).

### 5. Scoring, analyses and figures

```bash
python -m retrain.scoring.score main                 # every model of Table 1 and Figs. 1 and 3 -> results/main
python -m retrain.scoring.score selection            # Finch on the Fig. 4b control sets and the dose sets
python -m retrain.scoring.score explabels_training    # models trained on HQ with experimental labels only
python -m retrain.scoring.pr_curves                  # precision-recall curves
python -m retrain.scoring.test_labels                # scoring against experimental-only test labels
python -m retrain.scoring.no_hit                     # test proteins without training relatives (DIAMOND)
python -m retrain.scoring.sizes                      # sizes of every training set
python -m figures.make_all                           # every figure and table -> <data root>/figures
```

Scoring the main comparison takes about half a minute on one H100 (a few seconds per run on 32 CPU cores). The figure
scripts need only the result files; `--tikz` draws Fig. 2's PDF with true colour gradients through TikZ, as in the paper
(needs `tectonic` and the Arial font).

### External inputs

| Input | Source | Place under the data root | Checksum (SHA-256) |
|---|---|---|---|
| ArchipelaGO dataset release | Hugging Face, `jstempel/ArchipelaGO` (to be released) | `release/` | |
| DeepGO-SE training data | https://bio2vec.net/data/deepgo/deepgo2/training-data.tar.gz | `external/deepgose/` | `80cc5446ef2d32cfccd19ff6a634cb6c2946370caf0ec7b9fd1bfb3c2dc07a44` |
| GO-Plus (GO release 2026-07-26) | https://release.geneontology.org/2026-08-05/ontology/extensions/go-plus.owl | `external/deepgose/go-plus.owl` | `52f51e70a94667f887f89205daeb442abe5f6f8b6cd8dda51ff5ad7fe7524bf5` |
| TransFew release `TFewData.tar.gz` | https://doi.org/10.5281/zenodo.22139108 (21.0 GB) | `external/transfew/` (extracted) | `2f9007f96b1fdd28cb12c35e452f4415eaa69e8e6b2ae4e5026d100de672bf40` |
| DPFunc data archive (December 2025 version) | https://drive.google.com/file/d/1qrxbkk450GJzhVfqnAN9Ms798owrq96n | `external/dpfunc/dpfunc_published_dataset.zip` | `3fa11fc181b5db435c0ebe1ac6005ed38c1d0073283130b4b5767d163c1d56e2` |
| DPFunc InterPro vocabulary `inter_idx.pkl` | https://github.com/CSUBioGroup/DPFunc (commit `c5faf4df`), `data/inter_idx.pkl` | `external/dpfunc/` | |
| Language models | ESM2-15B, ESM2-3B, ESM-1b (fair-esm / Hugging Face), ProtT5-XL-U50 | downloaded on first use | |
| Programs | DIAMOND 2.2.4, InterProScan 5.78-109.0, ColabFold 1.6.3 with MMseqs2 release 18 and the ColabFold databases (UniRef30 2302, colabfold_envdb_202108), Java 11 with Groovy 2.4 | `PATH` | |

The scoring bundle written from the released benchmark should have the SHA-256 checksums that the paper's results record
(`inputs.json` of each results folder): `go-basic.obo` `b08d45b268b8c24ccb2513dbbbc7d4df9f6521c099b413f79eb31e06e0fa3bcc`,
`vocabulary.txt` `1fc60a7e2b8e7dd816baa9c804b4d48f61fc8efb73c86c077cd5652973018df4`.

## Layout

```
retrain/        common.py (datasets, ontology, fits, the homology screen), paths.py (the data root), and one package per
                stage: datasets, finch (training recipe, training, combination), baselines, deepgose, transfew, dpfunc,
                structures, scoring
scripts/        rank_coverage.py (the coverage ranking), embed.py (language-model embeddings), export_package.py (Finch's
                model package)
figures/        the paper's figures and tables (python -m figures.make_all)
configs/        finch.json: Finch's readouts, recipes, weights and seeds
data/           README.md: the data root's layout and file formats
examples/slurm/ an example job script
tests/          CPU tests: python -m pytest
```

## License

The code is released under the MIT License (`LICENSE`). Code adapted from DeepGO-SE (BSD 3-Clause), TransFew (MIT),
DPFunc (MIT) and fair-esm (MIT) keeps its license next to it; see `THIRD_PARTY_NOTICES.md`. The ArchipelaGO dataset is
released separately under CC BY-NC 4.0.
