"""DeepGO-SE pipeline logic that runs on a CPU: original-data parsing, early stopping, the 6-of-10 choice and the
ensemble's propagation (is_a and relationships, each world propagated before averaging)."""
import json

import numpy as np
import pandas as pd
import pytest

from retrain.common import Ontology
from retrain.deepgose import original, train

OBO = """format-version: 1.2

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0000002
name: a
namespace: molecular_function
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0000003
name: b
namespace: molecular_function
alt_id: GO:0000009
is_a: GO:0000002 ! a

[Term]
id: GO:0000004
name: c
namespace: molecular_function
is_a: GO:0003674 ! molecular_function
relationship: part_of GO:0000002 ! a

[Term]
id: GO:0000005
name: old
namespace: molecular_function
is_obsolete: true
replaced_by: GO:0000004

[Typedef]
id: part_of
name: part of
"""
VOCABULARY = ['GO:0000002', 'GO:0000003', 'GO:0000004']


@pytest.fixture
def root(tmp_path):
    folder = tmp_path / 'datasets' / 'ontology'
    folder.mkdir(parents=True)
    (folder / 'go-basic.obo').write_text(OBO)
    (folder / 'vocabulary.json').write_text(json.dumps(VOCABULARY))
    return tmp_path


def test_original_records_and_labels(root, tmp_path):
    source = tmp_path / 'mf'
    source.mkdir()
    for split, rows in (('train', [('P1', 'MKV', ['GO:0000009|IDA', 'GO:0005515|IEA'], ['GO:0000009'])]),
                        ('valid', [('P2', 'MKV', ['GO:0000005|EXP'], ['GO:0000005'])]),
                        ('test', [('P3', 'MAA', ['GO:0005515|IEA'], [])])):
        pd.DataFrame([dict(proteins=p, accessions=f'{p}_ACC;', orgs='9606', sequences=s, annotations=a, exp_annotations=e)
                      for p, s, a, e in rows]).to_pickle(source / f'{split}_data.pkl')
    proteins = original.merge(original.records(source))
    assert len(proteins) == 2  # P1 and P2 share a sequence
    merged = next(p for p in proteins.values() if p['sequence'] == 'MKV')
    assert merged['records'] == ['train:P1', 'valid:P2'] and merged['accessions'] == ['P1_ACC', 'P2_ACC']
    ontology = Ontology(root)
    # the alt_id maps to its primary term and the obsolete term to its single replacement
    assert original.mf_assertions(ontology, merged['codes']) == [('GO:0000005', 'GO:0000004', {'EXP'}),
                                                                ('GO:0000009', 'GO:0000003', {'IDA'})]
    other = next(p for p in proteins.values() if p['sequence'] == 'MAA')
    assert original.mf_assertions(ontology, other['codes']) == []  # no experimental annotation: left out


def test_should_stop():
    # best at epoch 4, patience 20: epochs 4..23 continue, epoch 24 (20 epochs without improvement) stops
    assert not any(train.should_stop(e, 4) for e in range(4, 24)) and train.should_stop(24, 4)
    assert train.should_stop(127, 120) and not train.should_stop(126, 120)  # the 128-epoch cap


def test_select_worlds():
    records = [dict(world=k, best_dev_bce=b) for k, b in enumerate([.5, .1, .3, .3, .9, .2, .8, .05, .7, .6])]
    assert train.select_worlds(records) == [7, 1, 5, 2, 3, 0]  # ties broken by world index


def test_ensemble_propagates_each_world_then_averages(root):
    rng = np.random.default_rng(0)
    scores, bces = {}, rng.permutation(10) / 100
    for k in range(10):
        folder = train.world_folder(root, 'hq', k)
        folder.mkdir(parents=True)
        for split, n in (('dev', 4), ('test', 5)):
            scores[k, split] = rng.random((n, 3), dtype=np.float32)
            np.save(folder / f'{split}_scores.npy', scores[k, split])
        (folder / 'record.json').write_text(json.dumps(dict(world=k, best_dev_bce=float(bces[k]), selected_epoch=1,
                                                            epochs_run=22)))
    train.ensemble(root, 'hq')
    seed0 = root / 'runs' / 'deepgose' / 'hq' / 'seed0'
    chosen = [int(k) for k in np.argsort(bces)[:6]]
    assert json.loads((seed0 / 'record.json').read_text())['ensemble_worlds'] == chosen

    def official(s):  # GO:0000002 takes the max of itself, its is_a child GO:0000003 and its part_of child GO:0000004
        out = s.copy()
        out[:, 0] = s.max(axis=1)
        return out
    for split in ('dev', 'test'):
        expected = np.mean([official(scores[k, split]).astype(np.float64) for k in chosen], axis=0).astype(np.float32)
        np.testing.assert_allclose(np.load(seed0 / f'{split}_scores.npy'), expected, rtol=0, atol=1e-7)
