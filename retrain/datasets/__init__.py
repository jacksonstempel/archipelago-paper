"""The training sets and the benchmark as dataset folders under <data root>/datasets/ (format in retrain/common.py).

release      the released dataset (Parquet) -> ontology, hq, archipelago, development, test
hq_subsets   goa (HQ proteins with GOA experimental or curator support, their GOA labels) and hq_explabels (HQ proteins
             with only the labels that make them HQ), both from datasets/hq
selected     sets built from the computational candidates: archipelago (check), archipelago_300k, archipelago_600k,
             pool_random, filter_random, coverage_unfiltered
random_draw  hq_random: HQ plus 150,000 proteins drawn at random from all computationally annotated TrEMBL proteins
candidates   datasets/candidates, the coverage ranking's candidates, from their release files
"""
