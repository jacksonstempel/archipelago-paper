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
pip install -e .[finch,embed,figures]  # Finch's training (GPU), embeddings, figures; add deepgose, transfew or dpfunc
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

The code here is the code that produced the paper's numbers, with the original runs' paths, logging and bookkeeping
removed. Three stages were rerun with it and reproduced the paper bit for bit on H100 GPUs (PyTorch 2.10.0 with CUDA 12.6,
Triton 3.6.0): Finch's training (`retrain/finch/engine.py`; run 9301 on ArchipelaGO, from the released data: both
checkpoints and every development and test score), the coverage ranking (both rankings, from the paper's input vectors)
and the embeddings (whole work units of ESM2-15B and ProtT5). The published methods are not bit-reproducible (GPU
nondeterminism), and where inputs that are not released change a result, the stage says so: the structure bank (section
4), the no-hit subsets (section 5) and the sets that need files available on request (data/README.md).

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
python scripts/rank_coverage.py --picks 600000                 # selection/coverage_ranking.tsv (four GPUs, minutes)
python scripts/rank_coverage.py --picks 500000 --no-filter     # selection/coverage_ranking_unfiltered.tsv (Fig. 3a control)
python -m retrain.datasets.selected check-archipelago          # the ranking's first 150,000 = the released selected proteins
python -m retrain.datasets.selected dose                       # archipelago_300k, archipelago_600k
python -m retrain.datasets.selected controls                   # pool_random, filter_random, coverage_unfiltered (DIAMOND)
python -m retrain.datasets.random_draw                         # hq_random, from its released protein list
```

The ranking is the paper's code (`archipelago.coverage`): run on the paper's candidate and HQ vectors on four H100 GPUs, it
reproduced both of the paper's rankings exactly. The paper's ranking read HQ vectors from an earlier extraction than the
released training embeddings (same encoder and windows; it handled X and U differently, and every row differs in the
last digits), so a ranking computed from the released embeddings can differ from the paper's among near-tied
candidates; the dataset builders use the paper's rankings instead when they are present (data/README.md).

### 3. Finch and the baselines

```bash
python scripts/embed.py --set archipelago development test --encoder esm2_15b_l48   # together, as the paper (GPU)
python scripts/embed.py --set archipelago development test --encoder prott5_l21     # then each further training set
python -m retrain.finch.train --training-set archipelago --seed 9301    # each training set, seeds 9301-9303 (GPU)
python -m retrain.finch.combine --training-set archipelago               # the 0.725 / 0.275 weighted mean
python -m retrain.baselines.naive --training-sets goa hq archipelago
python -m retrain.baselines.diamond_knn --training-sets goa hq archipelago        # DIAMOND
python -m retrain.baselines.embedding_knn --training-sets goa hq archipelago
python scripts/export_package.py --out finch_package           # the model package Finch.from_pretrained loads
```

Finch's recipe is `configs/finch.json`: an ESM2-15B readout (hidden width 8,192, dropout 0.8, ASAM rho 8) and a ProtT5
readout (width 4,096, dropout 0.6, rho 4), each trained with binary cross-entropy, AdamW without weight decay, 1,000
warm-up updates and a constant learning rate of 1e-5, ASAM eta 0.01 and batches of 256, with the best development
checkpoint kept (checks every 1,000 updates, patience 20); their raw scores are averaged 0.725 / 0.275. The package is
one ArchipelaGO run (seed 9301, `release_seeds`), with the threshold of its development micro-Fmax as the default cut.

The recipe ablation (Fig. 4a) builds Finch up on ArchipelaGO: a linear readout on ESM2-3B, the same on ESM2-15B, a narrow
hidden layer (width 1,024, dropout 0.3), the full hidden layer (width 8,192, dropout 0.8), ASAM (Finch's ESM2-15B readout),
and Finch. A control averages two ESM2-15B networks trained with different seeds; two comparisons train the full hidden
layer on ESM2-3B and on ESM-1b. Each model is listed under `ladder` in
`configs/finch.json` with the learning rate and dropout that the paper chose for it on the development set (the search
itself is not part of this repository). Finch with its ProtT5 readout on ProtT5's last layer (24) instead of layer 21,
same weights, is listed under `variants`:

```bash
python scripts/embed.py --set archipelago development test --encoder esm2_3b_l36   # (GPU)
python scripts/embed.py --set archipelago development test --encoder prott5_l24    # the same sets, for the variant
python -m retrain.finch.train --training-set archipelago --seed 9301 \
    --readouts esm2_3b_l36_linear esm2_15b_l48_linear esm2_15b_l48_narrow esm2_15b_l48_adamw esm2_15b_l48_second \
    esm2_3b_l36_adamw prott5_l24                                    # seeds 9301-9303
python -m retrain.finch.combine --training-set archipelago --model ladder_linear_esm2_3b       # each model under 'ladder'
python -m retrain.finch.combine --training-set archipelago --model finch_prott5_l24            # and under 'variants'
```

Embedding: ESM2-15B about 77 GPU-hours for the 9.9 million candidates, a few GPU-hours for the training sets; ProtT5 and
ESM2-3B are faster. `scripts/embed.py` is the paper's extraction (a fused FP16 ESM2 encoder, the paper's work units and
batches); embedding archipelago, development and test together reproduced the paper's vectors bit for bit (whole work
units of ESM2-15B and ProtT5 checked on an H100), while other groupings change the last digits (cosine 0.99997 or
closer). The environments it ran in are in its docstring. Finch needs a CUDA GPU and Triton (extra `finch`): 3 to 25 minutes per readout on one H100 (on
ArchipelaGO about 11 minutes for the ESM2-15B readout and 9 for the ProtT5 readout; about 3 GPU-hours for every training
set and seed of the main comparison). Baselines: minutes.

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
python -m retrain.structures.afdb requests && python -m retrain.structures.afdb fetch --chosen   # the paper's choice
# without inputs/structures/chosen.tsv: python -m retrain.structures.afdb search && python -m retrain.structures.afdb fetch
python -m retrain.structures.colabfold msa --databases DIR && python -m retrain.structures.colabfold fold
python -m retrain.structures.assemble
python -m retrain.dpfunc.embed                      # ESM-1b layer 31 per residue (GPU, about 4 GPU-hours)
python -m retrain.dpfunc.domains scan --interproscan DIR && python -m retrain.dpfunc.domains merge   # about 100 CPU-hours
python -m retrain.dpfunc.inputs graphs              # then: inputs bank <original|archipelago|evaluation>, inputs step <set>
python -m retrain.dpfunc.train --training-set <training set> --seed <seed>   # seeds 0-2; 23 GPU-hours for all 12 fits
```

Structures: an AlphaFold Database model of the exact sequence where one exists (84% of the paper's proteins),
otherwise an AlphaFold2 structure predicted with ColabFold (alphafold2_ptm, model 1, seed 0, 3 recycles, no templates).
`afdb fetch --chosen` downloads exactly the paper's AlphaFold Database entries from the paper's choice list
(data/README.md); a free search chooses differently for some proteins (the paper's bank preferred models already on
disk and kept 102 version-4 fragments of long human proteins). Most of the paper's 38,787 ColabFold structures were made
earlier in the project by other code, with faster kernels (a modified ColabFold) and mostly from MSAs of MMseqs2's GPU
search, neither of which this repository runs; a bank rebuilt here therefore differs from the paper's for many of those
proteins (`retrain/structures/colabfold.py` gives the counts).

**Finch on each method's original data.** Finch is trained on the three original sets with its recipe and weights
unchanged, as in section 3, once the sets exist:

```bash
python scripts/embed.py --set original_deepgose --encoder esm2_15b_l48     # and prott5_l21; original_transfew, original_dpfunc
python -m retrain.finch.train --training-set original_deepgose --seed 9301  # each set, seeds 9301-9303 (GPU)
python -m retrain.finch.combine --training-set original_deepgose
```

**The same network on DPFunc's language model.** The wide hidden layer without ASAM, on each protein's mean of DPFunc's
ESM-1b residue features:

```bash
python -m retrain.dpfunc.mean_features              # inputs/embeddings/esm1b_l31 (archipelago, development, test)
python -m retrain.finch.train --training-set archipelago --seed 9301 --readouts esm1b_l31_adamw   # seeds 9301-9303 (GPU)
python -m retrain.finch.combine --training-set archipelago --model ladder_hidden_esm1b
```

### 5. Scoring, analyses and figures

```bash
python -m retrain.scoring.score main                 # every model of Table 1 and Figs. 2 and 4b -> results/main
python -m retrain.scoring.score selection            # Finch on the Fig. 3a control sets and the dose sets
python -m retrain.scoring.score ladder               # the recipe ablation of Fig. 4a, its variants, Finch, published methods
python -m retrain.scoring.score explabels_training    # models trained on HQ with experimental labels only
python -m retrain.scoring.pr_curves                  # precision-recall curves on every shared training set (SI)
python -m retrain.scoring.distances                  # distances of the added proteins to HQ (Fig. 3b; GPU)
python -m retrain.scoring.test_labels experimental   # scoring against experimental-only test labels
python -m retrain.scoring.test_labels all_goa        # ... with every GOA assertion the label review considered
python -m retrain.scoring.test_labels masked         # ... with the reviewed GOA assertions masked
python -m retrain.scoring.term_strata                # the main gains on common and rare terms (GPU)
python -m retrain.scoring.shared_terms               # the data comparisons on the terms each original set covers (GPU)
python -m retrain.scoring.no_hit                     # test proteins without training relatives (Fig. 3c; the paper's subsets)
python -m retrain.scoring.sizes                      # sizes of every training set
python -m retrain.scoring.parameters                 # trainable parameters of the models on ArchipelaGO (CPU)
python -m retrain.scoring.training_cost              # GPU-hours of the training runs, from the fits' records
python -m figures.make_all                           # every figure and table -> <data root>/figures
```

Scoring the main comparison takes about half a minute on one H100 (a few seconds per run on 32 CPU cores). The figure
scripts need only the result files; `--tikz` draws Fig. 1's PDF with true colour gradients through TikZ, as in the paper
(needs `tectonic` and the Arial font).

### 6. The published methods on their own benchmarks (Supplementary Note 14)

**Replication.** Each published method is trained with this repository's code on its authors' own training split and scored
on their test set with their own evaluation code (`retrain.replicate`). A method's vocabulary and GO release are the
authors', so each method gets its own data root, here `rep/<method>`, with the authors' release under its `external/`
(the files of the table below; a link to the main data root's `external/` will do). The authors' splits share identical
sequences between training and test (TransFew 34, DPFunc 53), so their trainers run with `--allow-shared-proteins`; DPFunc's
77 vocabulary terms that are obsolete or merged in its GO release are listed in `datasets/ontology/not_live.json`.

```bash
python -m retrain.replicate.build deepgose --data rep/deepgose      # likewise transfew, dpfunc (CPU, minutes)
# DeepGO-SE: ESM2-3B inputs, ten worlds, the six-world ensemble, the released models, scoring
python -m retrain.deepgose.features original_deepgose development test --data rep/deepgose
python -m retrain.deepgose.train world original <k> --data rep/deepgose      # k = 0..9, one GPU each
python -m retrain.deepgose.train ensemble original --data rep/deepgose
python -m retrain.replicate.released deepgose --data rep/deepgose          # external/deepgose/data.tar.gz (GPU)
python -m retrain.replicate.score deepgose --testset --data rep/deepgose   # -> replicate/results, replicate/testset
# TransFew: ESM2-15B inputs, label GCN, seeds 0-2, scoring (its released prediction files are scored too)
python -m retrain.transfew.features embed --sets original_transfew development test --data rep/transfew
python -m retrain.transfew.features merge --sets original_transfew development test --data rep/transfew
python -m retrain.transfew.label_gcn --data rep/transfew
python -m retrain.transfew.train original <seed> --allow-shared-proteins --data rep/transfew
python -m retrain.replicate.score transfew --testset --data rep/transfew
# DPFunc: structures, ESM-1b, InterProScan, graphs and banks as in section 4 with --datasets original_dpfunc development
# test (and the banks original and evaluation), then seeds 0-2, the released checkpoints, scoring
python -m retrain.dpfunc.train --training-set original --seed <seed> --data rep/dpfunc
python -m retrain.replicate.released dpfunc --data rep/dpfunc              # external/dpfunc/save_models.tar.gz (GPU)
python -m retrain.replicate.score dpfunc --testset --data rep/dpfunc
```

The splits: DeepGO-SE 38,533 / 1,901 / 2,845 records (training, validation, test), TransFew 70,773 / 7,864 / 1,112 (Test_all;
Test_novel is 858 of them), DPFunc 31,463 / 682 / 401. Scoring follows each paper's protocol (docstring of
`retrain.replicate.score`): DeepGO-SE's `compute_metrics`, TransFew's CAFA-evaluator 1.2.1 call (`pip install
cafaeval==1.2.1`), DPFunc's `new_compute_performance_deepgoplus`; the authors' released DeepGO-SE models and TransFew
predictions reproduce their papers' numbers this way. The paper reports each method's Fmax (and Smin) by its own protocol
and, in place of the papers' protein-centric trapezoid AUPR, our micro-AUPRC on the matrices that protocol scores
(`retrain.replicate.micro`); `--paper-aupr` adds the papers' AUPR, the check against their tables. Costs on H100s: DeepGO-SE about 2 GPU-hours for the ten worlds,
TransFew 0.2 per seed, DPFunc about 0.8 per seed.

**Finch on the published methods' benchmarks.** Finch is trained on each method's training split and on ArchipelaGO
without the method's test proteins and their relatives (the paper's homology screen), with ArchipelaGO's
labels mapped onto the method's GO release and propagated over the relations the method's own training labels use
(`retrain.outside`). For the paper's runs, the proteins ArchipelaGO left out only because of the paper's own benchmark (its
2,500 development and 5,000 test proteins and their 4,063 HQ relatives, labelled as HQ training proteins) were added back
before the screen (`--addback`). The method's validation set is the development set of both arms (checkpoints and the F1
threshold); the test set is scored with the paper's own code. Both training sets can share sequences with the method's
validation set (and the authors' splits with its test set), hence `--allow-shared-proteins`.

```bash
python -m retrain.outside.build deepgose --data out/deepgose --replication rep/deepgose --archipelago <data root> [--addback NAME] [--screen-original]
python scripts/embed.py --set original_deepgose --encoder esm2_15b_l48 --data out/deepgose   # and development, test; prott5_l21
python -m retrain.finch.train --training-set original_deepgose --seed 9301 --allow-shared-proteins --data out/deepgose
python -m retrain.finch.train --training-set archipelago_screened --seed 9301 --allow-shared-proteins --data out/deepgose   # seeds 9301-9303 (GPU)
python -m retrain.finch.combine --training-set original_deepgose --data out/deepgose        # and archipelago_screened
python -m retrain.outside.score deepgose --data out/deepgose --replication rep/deepgose     # -> outside/results
```

`retrain.outside.build` needs the main data root's `datasets/archipelago` and its Finch embeddings
(`inputs/embeddings/{esm2_15b_l48,prott5_l21}/archipelago`, e.g. from the release as in data/README.md). The results are
reported for all test proteins and for those without a relative in either training set, and for the full vocabulary and
the terms with training positives in both arms. `--screen-original` adds a third training set, the method's split without
the test proteins' relatives under the same screen (`original_<method>_screened`; the paper: TransFew and DPFunc), trained
like `original_<method>` and scored with `retrain.outside.score --screened-original`; `compare` and `masked` include it
when its Finch runs exist.

One table per benchmark with every model (both Finch arms, the method's own retrain, the authors' model) and paired
bootstrap intervals, and the check of Supplementary Tables 16-17 with the test proteins' current GOA annotations masked
(our evaluator, micro-Fmax and micro-AUPRC; the bootstrap runs on a GPU):

```bash
python -m retrain.outside.compare deepgose --data out/deepgose --replication rep/deepgose   # -> outside/compare
mkdir -p <root>/results/outside_compare/deepgose                                            # figures/si.py's Table 18:
cp out/deepgose/outside/compare/relatives.tsv <root>/results/outside_compare/deepgose/      # likewise transfew, dpfunc
M="--data <root>/results/outside_masked --outside out --replication rep"                   # out/<method>, rep/<method>
python -m retrain.outside.masked fetch $M                                   # UniProt sec_ac.txt; DPFunc's entry names
python -m retrain.outside.masked scan $M --goa <GOA folder> --swissprot uniprot_sprot.fasta.gz   # goa_uniprot_all.{gaf,gpi}.gz
python -m retrain.outside.masked prepare deepgose $M                        # and transfew, dpfunc (CPU)
python -m retrain.outside.masked bootstrap deepgose $M                      # (GPU) -> Tables 16-17, Figs. 3-4
python -m retrain.outside.masked pooling deepgose $M                        # (GPU) -> Table 19
```

### External inputs

| Input | Source | Place under the data root | Checksum (SHA-256) |
|---|---|---|---|
| ArchipelaGO dataset release | Hugging Face, `jstempel/ArchipelaGO`, tag `v1.0` (to be released): `hf download jstempel/ArchipelaGO --repo-type dataset --revision v1.0 --local-dir <root>/release` | `release/` | `SHA256SUMS` in the release |
| DeepGO-SE training data | https://bio2vec.net/data/deepgo/deepgo2/training-data.tar.gz | `external/deepgose/` | `80cc5446ef2d32cfccd19ff6a634cb6c2946370caf0ec7b9fd1bfb3c2dc07a44` |
| GO-Plus (GO release 2026-07-26) | https://release.geneontology.org/2026-08-05/ontology/extensions/go-plus.owl | `external/deepgose/go-plus.owl` | `52f51e70a94667f887f89205daeb442abe5f6f8b6cd8dda51ff5ad7fe7524bf5` |
| TransFew release `TFewData.tar.gz` | https://doi.org/10.5281/zenodo.22139108 (21.0 GB) | `external/transfew/` (extracted) | `2f9007f96b1fdd28cb12c35e452f4415eaa69e8e6b2ae4e5026d100de672bf40` |
| DPFunc data archive (December 2025 version) | https://drive.google.com/file/d/1qrxbkk450GJzhVfqnAN9Ms798owrq96n | `external/dpfunc/dpfunc_published_dataset.zip` | `3fa11fc181b5db435c0ebe1ac6005ed38c1d0073283130b4b5767d163c1d56e2` |
| DPFunc InterPro vocabulary `inter_idx.pkl` | https://github.com/CSUBioGroup/DPFunc (commit `c5faf4df`), `data/inter_idx.pkl` | `external/dpfunc/` | |
| DeepGO-SE code (replication: its metrics) | https://github.com/bio-ontology-research-group/deepgo2, commit `6c77c98` | `external/deepgose/deepgo2/` | |
| DeepGO-SE trained models (replication) | https://bio2vec.net/data/deepgo/deepgo2/data.tar.gz (15.0 GB) | `external/deepgose/data.tar.gz` | |
| DPFunc code and label list (replication) | https://github.com/CSUBioGroup/DPFunc, commit `c5faf4df` | `external/dpfunc/DPFunc/` | |
| DPFunc trained models (replication) | https://drive.google.com/file/d/1V0VTFTiB29ilbAIOZn0okBQWPlbOI3wN (`save_models.tar.gz`, 14.3 GB) | `external/dpfunc/save_models.tar.gz` | |
| GO release 2022-12-04 (replication: DPFunc's labels) | https://release.geneontology.org/2022-12-04/ontology/go.obo | `external/go/2022-12-04/go.obo` | |
| Language models | ESM2-15B, ESM2-3B, ESM-1b (fair-esm / Hugging Face), ProtT5-XL-U50 | downloaded on first use | |
| Programs | DIAMOND 2.2.4, InterProScan 5.78-109.0, ColabFold 1.6.3 with MMseqs2 release 18 and the ColabFold databases (UniRef30 2302, colabfold_envdb_202108), Java 11 with Groovy 2.4 | `PATH` | |

The scoring bundle written from the released benchmark should have the SHA-256 checksums that the paper's results record
(`inputs.json` of each results folder): `go-basic.obo` `b08d45b268b8c24ccb2513dbbbc7d4df9f6521c099b413f79eb31e06e0fa3bcc`,
`vocabulary.txt` `1fc60a7e2b8e7dd816baa9c804b4d48f61fc8efb73c86c077cd5652973018df4`.

## Layout

```
retrain/        common.py (datasets, ontology, fits, the homology screen), paths.py (the data root), and one package per
                stage: datasets, finch (training recipe, training, combination), baselines, deepgose, transfew, dpfunc,
                structures, scoring, replicate (the methods on their own benchmarks), outside (Finch on them)
scripts/        rank_coverage.py (the coverage ranking), embed.py (language-model embeddings), export_package.py (Finch's
                model package)
figures/        the paper's figures and tables (python -m figures.make_all)
configs/        finch.json: Finch's readouts, recipes, weights and seeds, and the models of the recipe ablation
data/           README.md: the data root's layout and file formats
examples/slurm/ an example job script
tests/          CPU tests: python -m pytest
```

## License

The code is released under the MIT License (`LICENSE`). Code adapted from DeepGO-SE (BSD 3-Clause), TransFew (MIT),
DPFunc (MIT) and fair-esm (MIT) keeps its license next to it; see `THIRD_PARTY_NOTICES.md`. The ArchipelaGO dataset is
released separately under CC BY-NC 4.0.
