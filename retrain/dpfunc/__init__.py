"""DPFunc (Wang et al. 2025, Nature Communications) retrained with its published recipe.

Stages, in order (each module's --help gives details):
  original   datasets/original_dpfunc from the authors' release
  embed      ESM-1b layer-31 residue features (GPU)
  domains    InterProScan 5.78-109.0 domain counts (CPU)
  inputs     contact graphs from the structure bank (retrain.structures), the feature/edge banks, per-dataset rows
  train      one fit per training set and seed (GPU)
"""
