"""The paper's figures and tables, drawn from result files under the data root (python -m figures.make_all).

Files read, under <root>/results/ (each written by the retrain.scoring stage named):

  main/metrics.tsv, main/contrasts.tsv      score main          Table 1, Fig. 2, Fig. 4b's key, SI tables
  ladder/metrics.tsv, ladder/contrasts.tsv  score ladder        Fig. 4a, SI table of the recipe ablation
  pr_curves/curves.tsv                      pr_curves           Fig. 4b, SI precision-recall figure
  distances/densities.tsv, summary.tsv      distances           Fig. 3a
  selection/metrics.tsv, contrasts.tsv      score selection     Fig. 3b, SI tables of the selection rules and picks
  dataset_sizes.tsv                         sizes               Supplementary Table 2
  experimental_test_labels/                 test_labels         SI: experimental-only test labels
  explabels_training/                       score explabels_training   SI: HQ trained with experimental labels only
  no_hit/subsets.tsv, differences.tsv       no_hit              SI: test proteins without training relatives

Fig. 1 reads the counts in figures/sankey_counts.json; Extended Data Fig. 1 is a schematic.
DIAMOND-kNN is the method `diamond_knn` of the main run (no E-value cutoff, as reported).
"""
