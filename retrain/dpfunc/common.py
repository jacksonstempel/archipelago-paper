"""Constants and paths shared by the DPFunc stages.

<data root>/external/dpfunc/   the authors' inputs: dpfunc_published_dataset.zip (their benchmark data) and inter_idx.pkl
                               (the 22,369-entry InterPro vocabulary of the released model)
<data root>/inputs/dpfunc/     esm1b/ (residue features), domains/ (InterProScan counts), graphs/ (contact graphs),
                               banks/<bank>/ (features and edges in one file each), steps/<dataset>/ (bank rows, domains)
"""
from pathlib import Path

import numpy as np

N_INTERPRO = 22369
CONTACT_ANGSTROM = 10.0  # the paper (the released code's default is 12)
DIMENSION = 1280  # ESM-1b
TRAINING_SETS = ('original', 'goa', 'hq', 'archipelago', 'hq_explabels')  # hq_explabels: SI, experimental labels only
SEEDS = (0, 1, 2)
BANKS = {'original': ('original_dpfunc',), 'archipelago': ('archipelago',), 'evaluation': ('development', 'test')}
BANK_OF = {'original_dpfunc': 'original', 'archipelago': 'archipelago', 'goa': 'archipelago', 'hq': 'archipelago',
           'hq_explabels': 'archipelago', 'development': 'evaluation', 'test': 'evaluation'}  # goa, hq: subsets of archipelago


def dataset_of(training_set):
    """The dataset of a training set name ('original' is DPFunc's own data)."""
    return 'original_dpfunc' if training_set == 'original' else training_set


def bank_of(dataset):
    """The bank holding a dataset's proteins."""
    if dataset not in BANK_OF:
        raise ValueError(f'{dataset}: no DPFunc bank holds this dataset (known: {sorted(BANK_OF)})')
    return BANK_OF[dataset]


def external(root):
    return Path(root) / 'external' / 'dpfunc'


def inputs(root):
    return Path(root) / 'inputs' / 'dpfunc'


def esm_file(root, seq_id):
    return inputs(root) / 'esm1b' / seq_id[:2] / f'{seq_id}.npy'


def node_positions(length):
    """The residues that become graph nodes: all, or residues 1-1,400 above 2,700 (the structure bank's rule)."""
    return np.arange(1400 if length > 2700 else length, dtype=np.int64)
