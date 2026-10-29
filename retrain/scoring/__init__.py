"""Scoring of every fit with the paper's scorer (archipelago.scorer), and the analyses behind the Supplementary Information.

    python -m retrain.scoring.score RUN      the scorer on one set of fits (main, selection, ladder, explabels_training),
                                             then contrasts.tsv
    python -m retrain.scoring.pr_curves      precision-recall curves of Fig. 4b and the SI
    python -m retrain.scoring.distances      distances of the added proteins to HQ (Fig. 3a)
    python -m retrain.scoring.test_labels    scoring against experimental-only test labels
    python -m retrain.scoring.no_hit         scoring of test proteins without training relatives
    python -m retrain.scoring.sizes          sizes of the training sets and the benchmark

Every stage reads the data root (--data or $ARCHIPELAGO_DATA) and writes under <root>/results/<name>/; the figures
(figures/) read those folders.
"""
