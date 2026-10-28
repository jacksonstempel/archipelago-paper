"""TransFew (Boadu & Cheng, 2024) retrained with its published recipe on every training set.

Stages, in order (each `python -m retrain.transfew.<module> --help`):
  original   the authors' MF training data rebuilt from their release, screened against the benchmark (CPU, DIAMOND)
  features   ESM2-15B layer-48 inputs of every protein (GPU)
  label_gcn  the label graph over the vocabulary and its pretrained GCN (GPU, minutes)
  train      one fit per training set and seed (GPU)

model.py and label_gcn.py adapt code of the official repository (MIT; LICENSE-TransFew).
"""
