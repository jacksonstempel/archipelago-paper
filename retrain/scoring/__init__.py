"""Scoring of every fit with the paper's scorer (archipelago.scorer), and the analyses behind the Supplementary Information.

    python -m retrain.scoring.score RUN      the scorer on one set of fits (main, selection, ladder, explabels_training),
                                             then contrasts.tsv
    python -m retrain.scoring.pr_curves      precision-recall curves of Fig. 4b
    python -m retrain.scoring.distances      distances of the added proteins to HQ (Fig. 3a)
    python -m retrain.scoring.test_labels LABELS   scoring against other test labels: experimental-only, or every GOA
                                             assertion the label review considered (all_goa)
    python -m retrain.scoring.masked         scoring with the reviewed GOA assertions masked
    python -m retrain.scoring.no_hit         scoring of test proteins without training relatives
    python -m retrain.scoring.sizes          sizes of the training sets and the benchmark
    python -m retrain.scoring.parameters     trainable parameters of the models trained on ArchipelaGO
    python -m retrain.scoring.training_cost  GPU-hours of the training runs

Every stage reads the data root (--data or $ARCHIPELAGO_DATA) and writes under <root>/results/<name>/; the figures
(figures/) read those folders.
"""
