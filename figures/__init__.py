"""The paper's figures and tables, drawn from result files under the data root (python -m figures.make_all).

Files read, under <root>/results/ (each written by the retrain.scoring stage named):

  main/metrics.tsv, main/contrasts.tsv      score main          Tables 1 and 2, Figs. 2 and 4b, SI tables
  ladder/metrics.tsv, ladder/contrasts.tsv  score ladder        Fig. 4a, SI tables of the recipe ablation and paired
                                                                comparisons
  parameters.json                           parameters          SI table of the recipe ablation
  pr_curves/curves.tsv                      pr_curves           SI: precision-recall curves
  selection/metrics.tsv, contrasts.tsv      score selection     Fig. 3a, SI tables of the selection rules and picks
  distances/densities.tsv, summary.tsv      distances           Fig. 3b
  dataset_sizes.tsv                         sizes               Supplementary Table 2
  experimental_test_labels/                 test_labels experimental   Table 2, SI: the test versions
  all_goa_test_labels/                      test_labels all_goa        Table 2
  masked_test_labels/                       test_labels masked         Table 2
  direct_masked_test_labels/contrasts.tsv   label_grid          Table 2
  label_grid/grid.tsv, versions.tsv         label_grid          SI: every model on every version of the test set
  explabels_training/                       score explabels_training   Fig. 4b, SI: HQ with experimental labels only
  no_hit/subsets.tsv, differences.tsv,      no_hit              Table 2, Fig. 3c, SI: test proteins without training
    E<cutoff>/{metrics.tsv,contrasts.tsv}                       relatives
  term_strata/gains.tsv                     term_strata         SI: the main gains on common and rare terms
  macro/contrasts.tsv                       macro               SI: the main comparisons in term-centric AUPRC
  training_cost.tsv                         training_cost       SI: training cost
  outside_masked/                           retrain.outside.masked     SI: the published methods' own benchmarks

HEADLINE=masked (environment variable; default reviewed) draws every result that depends on the main test labels from 
<root>/results_masked/ instead (same folders and files; a missing file is an error) and shades the Annotations masked row of
Table 2 instead of the Reviewed labels row. Everything is then written to <root>/figures_masked/, so that the two sets of figures share no file. The label versions' own folders (experimental, direct, masked, all GOA) and the
no-hit subset sizes are read from results/ either way.

Fig. 1 reads the counts in figures/sankey_counts.json; Extended Data Fig. 1 is a schematic.
DIAMOND-kNN is the method `diamond_knn` of the main run (no E-value cutoff, as reported).
"""
