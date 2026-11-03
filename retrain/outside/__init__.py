"""Finch on the published methods' own benchmarks (Supplementary Note 14): on each method's test set, scored with its paper's
own code, Finch trained on the method's own training split against Finch trained on ArchipelaGO without the method's test
and validation proteins and their relatives.

  1. build.py   a data root per method: the replication root's ontology, `original_<method>`, `development` (the authors'
                validation set) and `test`, and `archipelago_screened`; the screen's counts and the test proteins' relatives
  2. scripts/embed.py for Finch's two inputs on the method's three sets; retrain.finch.train on both training sets (the
                original one with --allow-shared-proteins), retrain.finch.combine
  3. score.py   both arms with the paper's own code (retrain.replicate.score's protocols), F1 at each run's validation
                threshold, all test proteins and those without a relative in either training set, all and shared terms
"""
