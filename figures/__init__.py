"""The paper's figures and tables, drawn from result files under the data root (python -m figures.make_all).

Files read, under <root>/results/ (each written by the retrain.scoring stage named):

  main/metrics.tsv, main/contrasts.tsv      score main          Figs. 1 and 3, Table 1, SI tables
  pr_curves/curves.tsv                      pr_curves           Fig. 1b, SI precision-recall figure
  selection/metrics.tsv, contrasts.tsv      score selection     Fig. 4b, SI tables of the selection rules and picks
  dataset_sizes.tsv                         sizes               Supplementary Table 2
  experimental_test_labels/                 test_labels         SI: experimental-only test labels
  explabels_training/                       score explabels_training   SI: HQ trained with experimental labels only
  no_hit/subsets.tsv, differences.tsv       no_hit              SI: test proteins without training relatives

Fig. 2 reads the counts in figures/sankey_counts.json; Fig. 4a and Extended Data Fig. 1 are schematics.
DIAMOND-kNN is the method `diamond_knn` of the main run (no E-value cutoff, as reported).
"""
