"""What retrain.scoring adds to the scorer (whose end-to-end test is tests/test_scorer.py): the scoring bundle exported
from the dataset folders, read back as archipelago.scorer reads it; contrasts.tsv = the difference of the scorer's
seed means with the paired interval of its replicates; and the masked scoring, which equals the scorer's metrics without a
mask and weighs proteins by their cluster draws in the bootstrap (CPU, tiny synthetic data)."""
import csv
import json

import numpy as np
import pytest

pytest.importorskip('scipy')
pytest.importorskip('torch')

import torch  # noqa: E402

from archipelago import evaluate  # noqa: E402
from archipelago.evaluate import paired_difference  # noqa: E402
from archipelago.io import read_labels, read_table, read_vocabulary  # noqa: E402
from retrain.common import Ontology, load, seq_id, write_dataset  # noqa: E402
from retrain.scoring import score  # noqa: E402
from retrain.scoring.masked import metrics, replicate_auprc  # noqa: E402

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


def masked_case():
    """Scores with ties and on the thresholds (multiples of 1/40), labels, information accretion."""
    rng = np.random.default_rng(2)
    scores = (rng.integers(0, 41, (30, 12)) / 40).astype(np.float32)
    labels = rng.random((30, 12)) < 0.3
    return rng, scores, labels, rng.random(12) * 3


def test_masked_metrics_without_a_mask_are_the_scorers():
    rng, s, y, ia = masked_case()
    as_tensor = torch.as_tensor
    found = metrics(as_tensor(s), as_tensor(y), as_tensor(np.ones_like(y)), as_tensor(ia))
    assert found['auprc'] == pytest.approx(evaluate.micro_auprc(s, y), abs=1e-12)
    assert found['fmax'] == pytest.approx(evaluate.micro_fmax(s, y)[0], abs=1e-12)
    assert found['smin'] == pytest.approx(evaluate.smin(s, y, ia)[0], abs=1e-12)
    keep = y | (rng.random(y.shape) < 0.6)  # masked pairs count as if they did not exist
    found = metrics(as_tensor(s), as_tensor(y), as_tensor(keep), as_tensor(ia))
    assert found['auprc'] == pytest.approx(evaluate.micro_auprc(s[keep][None], y[keep][None]), abs=1e-12)
    assert found['fmax'] == pytest.approx(evaluate.micro_fmax(s[keep][None], y[keep][None])[0], abs=1e-12)


def test_masked_bootstrap_weighs_proteins_by_their_cluster_draws():
    rng, s, y, _ = masked_case()
    cluster = np.arange(30) // 3  # ten clusters of three proteins
    draws = rng.integers(0, 3, (10, 5))
    draws[0] = 1  # every replicate holds positives
    keep = y | (rng.random(y.shape) < 0.6)
    found = replicate_auprc(torch.as_tensor(s), torch.as_tensor(y), torch.as_tensor(keep), torch.as_tensor(cluster),
                            torch.as_tensor(draws, dtype=torch.float64), chunk=2)
    for r in range(draws.shape[1]):
        rows = np.repeat(np.arange(30), draws[cluster, r])  # each protein as often as its cluster was drawn
        expected = evaluate.micro_auprc(s[rows][keep[rows]][None], y[rows][keep[rows]][None])
        assert found[r] == pytest.approx(expected, abs=1e-12)
