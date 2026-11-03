"""retrain.finch (train, reload, combine) and the naive and embedding-kNN baselines on a tiny fake data root. Training
needs a CUDA GPU and Triton (the paper's engine) and is skipped without them; the input checks and the combination run on
a CPU."""
import hashlib
import json
import sys

import numpy as np
import pytest
import torch

from archipelago.io import write_matrix
from retrain import common
from retrain.baselines import embedding_knn, naive
from retrain.finch import combine, config
from retrain.finch.train import Inputs, fit, recipe_of

try:
    import triton  # noqa: F401
    GPU = torch.cuda.is_available()
except ImportError:
    GPU = False
needs_gpu = pytest.mark.skipif(not GPU, reason='the Finch engine needs a CUDA GPU and Triton')

OBO = """[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

""" + ''.join(f"""[Term]
id: GO:000000{i}
name: t{i}
namespace: molecular_function
is_a: GO:0003674 ! molecular_function

""" for i in range(1, 5))
TERMS = [f'GO:000000{i}' for i in range(1, 5)]


def fake_root(root, sizes=dict(archipelago=60, development=12, test=8)):
    folder = root / 'datasets' / 'ontology'
    folder.mkdir(parents=True)
    (folder / 'go-basic.obo').write_text(OBO, encoding='utf-8')
    (folder / 'vocabulary.json').write_text(json.dumps(TERMS), encoding='utf-8')
    ontology = common.Ontology(root)
    rng = np.random.default_rng(0)
    for name, n in sizes.items():
        sequences = [f'M{name[:3].upper()}{i}'.replace('0', 'A').replace('1', 'C').replace('2', 'D')
                     .replace('3', 'E').replace('4', 'F').replace('5', 'G').replace('6', 'H').replace('7', 'I')
                     .replace('8', 'K').replace('9', 'L') for i in range(n)]
        ids = [hashlib.sha256(s.encode()).hexdigest() for s in sequences]
        terms = rng.integers(0, 4, n)
        rows = [(h, TERMS[t], 'IEA', 'GOA:test') for h, t in zip(ids, terms)]
        labels, _ = ontology.label_matrix([(h, g) for h, g, _, _ in rows], {h: i for i, h in enumerate(ids)})
        common.write_dataset(root, name, ids, [[] for _ in ids], sequences, rows, labels, ontology)
        for readout, dim in (('esm2_15b_l48', 6), ('prott5_l21', 4)):  # embeddings that carry the label, plus noise
            x = rng.random(size=(n, dim)).astype(np.float32)  # positive, like real embeddings' cosines
            x[np.arange(n), terms] += 3
            write_matrix(root / 'inputs' / 'embeddings' / readout / name, ids, x.astype(np.float16))


def small(settings, readout):
    """A readout's recipe at the scale of the fake data root."""
    return dict(recipe_of(settings, readout), width=8, batch=8, warmup=2, check_interval=5, patience=2)


@needs_gpu
def test_finch_train_and_combine(tmp_path, monkeypatch):
    fake_root(tmp_path)
    settings = config()
    for readout in settings['weights']:
        inputs = Inputs(tmp_path, 'archipelago', settings['readouts'][readout]['embedding'])
        fit(tmp_path, 'archipelago', 9301, readout, small(settings, readout), inputs)
    member = tmp_path / 'runs' / 'finch' / 'archipelago' / 'seed9301' / 'esm2_15b_l48'
    saved = torch.load(member / 'readout.pt', weights_only=True)
    assert saved['vocabulary'] == TERMS and saved['input_dim'] == 6
    assert np.load(member / 'test_scores.npy').shape == (8, 4)
    history = json.loads((member / 'history.json').read_text(encoding='utf-8'))
    assert history['reproduced_dev_max_abs'] <= 1e-5 and history['seed'] == 9301
    monkeypatch.setattr(sys, 'argv', ['x', '--training-set', 'archipelago', '--seeds', '9301', '--data', str(tmp_path)])
    combine.main()
    folder = member.parent
    expected = sum(np.float32(w) * np.load(folder / r / 'dev_scores.npy') for r, w in settings['weights'].items())
    assert np.array_equal(np.load(folder / 'dev_scores.npy'), expected.astype(np.float32))
    assert json.loads((folder / 'record.json').read_text(encoding='utf-8'))['weights'] == settings['weights']

    # the recipe ablation's control: a second ESM2-15B network trained with seed 9301 + 3, averaged with the first
    inputs = Inputs(tmp_path, 'archipelago', 'esm2_15b_l48')
    fit(tmp_path, 'archipelago', 9301, 'esm2_15b_l48_second', small(settings, 'esm2_15b_l48_second'), inputs, 3)
    assert torch.load(folder / 'esm2_15b_l48_second' / 'readout.pt', weights_only=True)['seed'] == 9304
    monkeypatch.setattr(sys, 'argv', ['x', '--training-set', 'archipelago', '--seeds', '9301', '--model',
                                      'ladder_two_esm2_15b', '--data', str(tmp_path)])
    combine.main()
    pair = tmp_path / 'runs' / 'ladder_two_esm2_15b' / 'archipelago' / 'seed9301'
    expected = sum(np.float32(0.5) * np.load(folder / r / 'dev_scores.npy') for r in ('esm2_15b_l48', 'esm2_15b_l48_second'))
    assert np.array_equal(np.load(pair / 'dev_scores.npy'), expected.astype(np.float32))


def fake_readout(root, readout, seed, rng, training_set='archipelago'):
    """A trained readout's files as retrain.finch.train writes them (random scores)."""
    folder = common.fit_folder(root, 'finch', training_set, 9301) / readout
    folder.mkdir(parents=True)
    torch.save(dict(state_dict={}), folder / 'readout.pt')
    for split, n in (('dev', 12), ('test', 8)):
        np.save(folder / f'{split}_scores.npy', rng.random((n, 4), dtype=np.float32))
    (folder / 'history.json').write_text(json.dumps(dict(readout=readout, training_set=training_set, seed=seed,
                                                         best_step=1000, dev_metrics={})), encoding='utf-8')
    return folder


def test_combine_is_the_weighted_sum_of_the_readouts(tmp_path, monkeypatch):
    """The paper's expression (float32 weight x float32 raw scores, summed in the listed order), and a readout of another
    seed is refused (CPU)."""
    rng = np.random.default_rng(0)
    settings = config()
    folders = {r: fake_readout(tmp_path, r, 9301, rng) for r in settings['weights']}
    monkeypatch.setattr(sys, 'argv', ['x', '--training-set', 'archipelago', '--seeds', '9301', '--data', str(tmp_path)])
    combine.main()
    for split in ('dev', 'test'):
        expected = np.float32(0) + sum(np.float32(w) * np.load(folders[r] / f'{split}_scores.npy')
                                       for r, w in settings['weights'].items())
        assert np.array_equal(np.load(folders['esm2_15b_l48'].parent / f'{split}_scores.npy'), expected)
    fake_readout(tmp_path, 'esm2_15b_l48_second', 9301, rng)  # trained with the run's seed, not seed + 3
    monkeypatch.setattr(sys, 'argv', ['x', '--training-set', 'archipelago', '--seeds', '9301', '--model',
                                      'ladder_two_esm2_15b', '--data', str(tmp_path)])
    with pytest.raises(ValueError, match='expected'):
        combine.main()


def test_finch_refuses_training_proteins_in_the_benchmark(tmp_path):
    fake_root(tmp_path, sizes=dict(development=12, test=8))
    dev = common.load(tmp_path, 'development')
    ontology = common.Ontology(tmp_path)
    common.write_dataset(tmp_path, 'leaky', dev.seq_ids, [[] for _ in dev.seq_ids], dev.sequences, [], dev.labels, ontology)
    with pytest.raises(ValueError, match='development or test'):  # before any GPU work
        Inputs(tmp_path, 'leaky', 'prott5_l21')


def test_naive_and_embedding_knn(tmp_path, monkeypatch):
    fake_root(tmp_path)
    for module in (naive, embedding_knn):
        monkeypatch.setattr(sys, 'argv', ['x', '--training-sets', 'archipelago', '--data', str(tmp_path)]
                            + (['--device', 'cpu'] if module is embedding_knn else []))
        module.main()
    folder = tmp_path / 'runs' / 'naive' / 'archipelago' / 'seed0'
    frequency = common.load(tmp_path, 'archipelago').labels.toarray().mean(0)
    assert np.allclose(np.load(folder / 'test_scores.npy'), frequency[None], atol=1e-6)
    record = json.loads((tmp_path / 'runs' / 'embedding_knn' / 'archipelago' / 'seed0' / 'record.json').read_text(encoding='utf-8'))
    assert record['k'] in (1, 3, 5, 10, 20, 50) and record['dev_auprc'] > 0.5  # the embeddings carry the labels
