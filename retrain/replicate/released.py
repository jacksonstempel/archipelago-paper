"""The authors' released models on their own test set, with the replication root's test inputs.

    python -m retrain.replicate.released deepgose --data DIR     (one GPU) -> replicate/released/test_scores.npy
    python -m retrain.replicate.released dpfunc --data DIR       (one GPU) -> replicate/released/test_scores_raw.npy

DeepGO-SE  external/deepgose/data.tar.gz (https://bio2vec.net/data/deepgo/deepgo2/data.tar.gz, 15,046,938,364 bytes): its
           MF models 0, 1, 2, 5, 6, 8 (predict.py's MF ensemble; deepgozero_esm_plus_<k>.th), each run on the test
           proteins' ESM2-3B features (retrain.deepgose.features: the official extract_esm computation), propagated over
           go.obo with every relationship (as train.py's propagate_annots) and averaged (evaluate_entailment.py -c avg).
DPFunc     external/dpfunc/save_models.tar.gz (the README's "Trained Model Weights", Google Drive
           1V0VTFTiB29ilbAIOZn0okBQWPlbOI3wN, 14,268,592,581 bytes): its three MF checkpoints
           (DPFunc_model_mf_<k>of3model.pt, DPFunc_main.py's three-slot rule), each run on the test proteins' inputs of this
           root (ESM-1b, 10 Angstrom contacts, InterProScan 5.78; the released weights may have been trained with ESM2-650M and
           12 Angstrom, as the released code uses), and the mean of their sigmoid scores (merge_result).
TransFew   its weights were part of TFewData; its released prediction files are scored directly (score.py).
The needed members are extracted once into replicate/released/; rows follow datasets/test. score.py scores them.
"""
import argparse
import fnmatch
import json
import tarfile
from pathlib import Path

import numpy as np

from retrain.common import add_argument, data_root, load, log

DEEPGOSE_MODELS = (0, 1, 2, 5, 6, 8)


def extract(archive, folder, wanted):
    """Members of a .tar.gz whose name matches one of the `wanted` patterns -> folder (flat); returns their paths."""
    folder.mkdir(parents=True, exist_ok=True)
    found = []
    with tarfile.open(archive, 'r|gz') as tar:
        for member in tar:
            name = Path(member.name).name
            if member.isfile() and any(fnmatch.fnmatch(member.name.lstrip('./'), w) for w in wanted):
                (folder / name).write_bytes(tar.extractfile(member).read())
                found.append(folder / name)
    return sorted(found)


def deepgose(root):
    import pandas as pd
    import torch
    from archipelago.io import read_matrix
    from retrain.deepgose.model import DeepGOModel, load_normal_forms
    from retrain.deepgose.train import BATCH, OfficialPropagation
    folder = root / 'replicate' / 'released'
    wanted = ['data/go.obo', 'data/go-plus.norm', 'data/mf/terms.pkl',
              *(f'data/mf/deepgozero_esm_plus_{k}.th' for k in DEEPGOSE_MODELS)]
    if not all((folder / Path(w).name).exists() for w in wanted):
        extract(root / 'external' / 'deepgose' / 'data.tar.gz', folder, wanted)
    terms = list(pd.read_pickle(folder / 'terms.pkl')['gos'].values.flatten())
    if terms != json.loads((root / 'datasets/ontology/vocabulary.json').read_text()):
        raise ValueError("the released models' terms.pkl differs from the root's vocabulary")
    forms = load_normal_forms(str(folder / 'go-plus.norm'), {t: i for i, t in enumerate(terms)})
    test = load(root, 'test', labels=False, sequences=False)
    _, x = read_matrix(root / 'inputs/deepgose/features/test', test.seq_ids, mmap=False)
    x = torch.from_numpy(np.asarray(x, dtype=np.float32)).cuda()
    torch.backends.cuda.matmul.allow_tf32 = False
    propagate = OfficialPropagation(folder / 'go.obo', terms)
    total = np.zeros((len(test.seq_ids), len(terms)), dtype=np.float64)
    for k in DEEPGOSE_MODELS:
        model = DeepGOModel(2560, len(terms), len(forms[5]), len(forms[4]), 'cuda:0').cuda()
        model.load_state_dict(torch.load(folder / f'deepgozero_esm_plus_{k}.th', map_location='cuda'))
        model.eval()
        with torch.inference_mode():
            raw = torch.cat([model(x[i:i + BATCH]) for i in range(0, len(x), BATCH)]).cpu().numpy()
        total += propagate(raw)
        log(f'released model {k} scored')
    np.save(folder / 'test_scores.npy', (total / len(DEEPGOSE_MODELS)).astype(np.float32))


def dpfunc(root):
    import torch
    from retrain.dpfunc.train import Split, evaluate, new_model
    folder = root / 'replicate' / 'released'
    checkpoints = sorted(folder.glob('*_mf_*of3model.pt'))
    if not checkpoints:
        checkpoints = extract(root / 'external' / 'dpfunc' / 'save_models.tar.gz', folder, ['*_mf_*of3model.pt'])
    if len(checkpoints) != 3:
        raise SystemExit(f'expected three MF checkpoints, found {[p.name for p in checkpoints]}')
    vocabulary = json.loads((root / 'datasets/ontology/vocabulary.json').read_text())
    test = Split(root, 'test', len(vocabulary), labels=False)
    parts = []
    for path in checkpoints:
        model = new_model(len(vocabulary), 'cuda')
        model.load_state_dict(torch.load(path, map_location='cuda')['model_state_dict'])
        parts.append(evaluate(model, test)[0])
        log(f'{path.name} scored')
    np.save(folder / 'test_scores_raw.npy', np.mean(np.stack(parts), axis=0).astype(np.float32))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('method', choices=['deepgose', 'dpfunc'])
    add_argument(parser)
    a = parser.parse_args()
    globals()[a.method](data_root(a.data))


if __name__ == '__main__':
    main()
