"""The plain network on DPFunc's language model: each protein's mean of DPFunc's ESM-1b residue features, as the input of
Finch's trainer (the readout esm1b_l31_adamw of configs/finch.json, reported in the Results and the Supplementary Information).

    python -m retrain.dpfunc.mean_features [--data DIR]

For ArchipelaGO, the development set and the test set: the mean over all of a protein's residues of inputs/dpfunc/esm1b
(retrain.dpfunc.embed), accumulated in float64, written to inputs/embeddings/esm1b_l31/<set>.{npy,ids.txt} (FP16, rows in
dataset order; the .ids.txt file is written last).
"""
import argparse

import numpy as np

from retrain.common import add_argument, data_root, load, log
from retrain.dpfunc.common import esm_file

SETS = ('archipelago', 'development', 'test')
DIMENSIONS = 1280


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argument(parser)
    root = data_root(parser.parse_args().data)
    folder = root / 'inputs' / 'embeddings' / 'esm1b_l31'
    folder.mkdir(parents=True, exist_ok=True)
    for name in SETS:
        ids = load(root, name, labels=False, sequences=False).seq_ids
        stored = np.lib.format.open_memmap(folder / f'{name}.npy', mode='w+', dtype=np.float16, shape=(len(ids), DIMENSIONS))
        for i, seq_id in enumerate(ids):
            stored[i] = np.load(esm_file(root, seq_id)).astype(np.float64).mean(0)
        stored.flush()
        (folder / f'{name}.ids.txt').write_text('\n'.join(ids) + '\n', encoding='utf-8')
        log(f'{name}: {len(ids):,} proteins')


if __name__ == '__main__':
    main()
