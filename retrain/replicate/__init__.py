"""The published methods on their own benchmarks (Supplementary Note 14): each method trained with this repository's code
on its authors' own training split and scored on their test set with their own evaluation code.

One data root per method (its vocabulary and GO release are the authors'), with the authors' release under external/:

  1. build.py      the authors' train / validation / test split as datasets/{ontology, original_<method>, development,
                   test}, and replicate/records.tsv (each authors' record and its protein)
  2. the method's own stages (README) compute its inputs and train it on `original`, with --allow-shared-proteins where the
     authors' split shares identical sequences between training and test (TransFew, DPFunc)
  3. score.py      each paper's own evaluation code on the test set; also writes the test set as scored (replicate/testset)
  4. released.py   the authors' released models, run on the same test inputs and scored the same way
"""
