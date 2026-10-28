"""What retrain.scoring adds to the scorer (whose end-to-end test is tests/test_scorer.py): the scoring bundle exported
from the dataset folders, read back as archipelago.scorer reads it, and contrasts.tsv = the difference of the scorer's
seed means with the paired interval of its replicates (CPU, tiny synthetic data root)."""
import csv
import json

import numpy as np
import pytest

pytest.importorskip('scipy')
pytest.importorskip('torch')

from archipelago.evaluate import paired_difference  # noqa: E402
from archipelago.io import read_labels, read_table, read_vocabulary  # noqa: E402
from retrain.common import Ontology, load, seq_id, write_dataset  # noqa: E402
from retrain.scoring import score  # noqa: E402

OBO = """format-version: 1.2

[Term]
id: GO:0003674
name: molecular_function

[Term]
id: GO:0000001
name: a
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0000002
name: b
is_a: GO:0000001 ! a

[Term]
id: GO:0000003
name: c
is_a: GO:0003674 ! molecular_function
"""
VOCABULARY = ['GO:0000001', 'GO:0000002', 'GO:0000003']


def tiny_root(root):
    rng = np.random.default_rng(0)
    (root / 'datasets' / 'ontology').mkdir(parents=True)
    (root / 'datasets' / 'ontology' / 'go-basic.obo').write_text(OBO)
    (root / 'datasets' / 'ontology' / 'vocabulary.json').write_text(json.dumps(VOCABULARY))
    ontology = Ontology(root)
    for name, n in (('hq', 12), ('development', 6), ('test', 10)):
        sequences = [''.join(rng.choice(list('ACDEFGHIKLMNPQRSTVWY'), 30)) for _ in range(n)]
        ids = [seq_id(s) for s in sequences]
        terms = [['GO:0000002', 'GO:0000003'][i % 2] for i in range(n)]
        labels, _ = ontology.label_matrix(list(zip(ids, terms)), {h: i for i, h in enumerate(ids)})
        columns = {'cluster': [str(i // 2) for i in range(n)]} if name != 'hq' else None
        write_dataset(root, name, ids, [[] for _ in ids], sequences, [(h, t, 'EXP', 'GOA:X') for h, t in zip(ids, terms)],
                      labels, ontology, columns=columns)


def test_bundle_reads_back_as_the_datasets(tmp_path):
    tiny_root(tmp_path)
    data = score.bundle(tmp_path)
    vocabulary = read_vocabulary(data / 'vocabulary.txt')
    assert vocabulary == VOCABULARY
    for split, table, labels in (('development', 'benchmark/development.tsv', 'benchmark/development_labels.tsv'),
                                 ('test', 'benchmark/test.tsv', 'benchmark/test_labels.tsv'),
                                 ('hq', 'archipelago/proteins.tsv', 'archipelago/labels.tsv')):
        ds = load(tmp_path, split, scoring=True)
        rows = read_table(data / table)
        assert [r['protein_id'] for r in rows] == ds.seq_ids
        if split == 'hq':
            assert {r['part'] for r in rows} == {'hq'}  # the scorer estimates information accretion from these
        else:
            assert [r['cluster'] for r in rows] == ds.columns['cluster']
        assert (read_labels(data / labels, ds.seq_ids, vocabulary) == ds.labels.toarray()).all()


def test_contrasts_are_paired_differences_of_seed_means(tmp_path):
    rng = np.random.default_rng(1)
    metrics = ['auprc', 'fmax', 'f1', 'smin', 'protein_fmax']
    models = [('finch', 'hq'), ('finch', 'archipelago'), ('naive', 'hq')]
    replicates = {f'{m}/{s}/mean': rng.random((200, len(metrics))) for m, s in models}
    np.savez(tmp_path / 'replicates.npz', metrics=np.array(metrics), **replicates)
    point = {(m, s, k): float(rng.random()) for m, s in models for k in metrics}
    with open(tmp_path / 'metrics.tsv', 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['method', 'training_set', 'seed', 'split', 'metric', 'value'])
        for (m, s, k), v in point.items():
            writer.writerow([m, s, 'mean', 'test', k, repr(v)])
            writer.writerow([m, s, 'mean', 'development', k, '0.5'])
    rows = {(r['metric'], r['method'], r['training_set'], r['minus_method'], r['minus_training_set']): r
            for r in score.contrasts(tmp_path)}
    assert len(rows) == 4 * (2 + 2)  # per metric: both orders of the finch pair and of the pair on hq
    for (m1, s1), (m2, s2) in (((('finch', 'archipelago'), ('finch', 'hq'))), (('naive', 'hq'), ('finch', 'hq'))):
        for k in ('auprc', 'smin'):
            j = metrics.index(k)
            expected, (low, high) = paired_difference((point[m1, s1, k], replicates[f'{m1}/{s1}/mean'][:, j]),
                                                      (point[m2, s2, k], replicates[f'{m2}/{s2}/mean'][:, j]))
            r = rows[k, m1, s1, m2, s2]
            assert r['difference'] == pytest.approx(expected)
            assert (r['ci_low'], r['ci_high']) == pytest.approx((low, high))
