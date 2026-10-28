"""DeepGO-SE retrained with its published recipe on every training set.

Stage order (each module's docstring has the details):
  1. original.py         datasets/original_deepgose from the authors' processed training data (CPU, DIAMOND)
  2. axioms/build.sh     GO-Plus axioms, EL profile, normalized (Java, Groovy)
  3. features.py         ESM2-3B layer-36 inputs of every dataset (one GPU)
  4. train.py world      ten worlds per training set (one GPU each)
  5. train.py ensemble   the six best worlds, propagated and averaged (CPU)

model.py, train.py's official_ancestors (deepgo/utils.py), makeEl.groovy and most of Normalizer.groovy come from
github.com/bio-ontology-research-group/deepgo2 (BSD 3-Clause; LICENSE-deepgo2).
"""
