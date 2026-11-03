"""What retrain.scoring adds to the scorer (whose end-to-end test is in the archipelago package): the scoring bundle exported
from the dataset folders, read back as archipelago.scorer reads it; contrasts.tsv = the difference of the scorer's
seed means with the paired interval of its replicates; the versions of the test labels, rebuilt from the annotations,
checked against the release and scored; and the bootstrap of term_strata, which weighs proteins by their cluster draws
(CPU, tiny synthetic data)."""
import csv
import gzip
import json
import sys

import numpy as np
import pytest

pytest.importorskip('scipy')
pytest.importorskip('torch')
pa = pytest.importorskip('pyarrow')
pq = pytest.importorskip('pyarrow.parquet')

import scipy.sparse as sp  # noqa: E402
import torch  # noqa: E402

from archipelago import evaluate  # noqa: E402
from archipelago.evaluate import paired_difference  # noqa: E402
from archipelago.io import read_labels, read_table, read_vocabulary  # noqa: E402
from retrain.common import Ontology, load, seq_id, write_dataset  # noqa: E402
from retrain.scoring import score, test_labels  # noqa: E402
from retrain.scoring.term_strata import replicate_auprc  # noqa: E402

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


def test_contrasts_set_finch_beside_each_method_on_its_original_data(tmp_path):
    rng = np.random.default_rng(3)
    metrics = ['auprc', 'fmax', 'f1', 'smin']
    models = [('finch', 'original_dpfunc'), ('dpfunc', 'original')]
    np.savez(tmp_path / 'replicates.npz', metrics=np.array(metrics),
             **{f'{m}/{s}/mean': rng.random((50, len(metrics))) for m, s in models})
    point = {('finch', 'original_dpfunc'): 0.45, ('dpfunc', 'original'): 0.40}
    with open(tmp_path / 'metrics.tsv', 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['method', 'training_set', 'seed', 'split', 'metric', 'value'])
        writer.writerows([m, s, 'mean', 'test', k, repr(v)] for (m, s), v in point.items() for k in metrics)
    rows = score.contrasts(tmp_path)
    assert len(rows) == len(score.CONTRAST_METRICS)  # no shared set, no shared model: only the cross pair
    r = next(r for r in rows if r['metric'] == 'auprc')
    assert (r['method'], r['training_set'], r['minus_method'], r['minus_training_set']) == ('finch', 'original_dpfunc', 'dpfunc', 'original')
    assert r['difference'] == pytest.approx(0.05)


def release_of(root):
    """release/ for a tiny data root, as Hugging Face holds it, and one fit. Test proteins: the first has only the review's
    accepted assertion besides a curator statement, so it has no experimental label; the second has an assertion that the
    review did not accept."""
    (root / 'release' / 'ontology').mkdir(parents=True)
    with gzip.open(root / 'release' / 'ontology' / 'go-basic.obo.gz', 'wt') as handle:
        handle.write(OBO)

    def write(name, rows):
        path = root / 'release' / f'{name}.parquet'
        path.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.Table.from_pylist(rows), path)

    write('go_terms/go_terms', [dict(go_id=t, column=j) for j, t in enumerate(VOCABULARY)])
    hq = load(root, 'hq', scoring=True)
    y = hq.labels.toarray()
    write('training/train', [dict(protein_id=h, accessions=[], part='hq', go_terms=[t for j, t in enumerate(VOCABULARY) if y[i, j]])
                             for i, h in enumerate(hq.seq_ids)])
    ontology = Ontology(root)
    fits = {}
    for split in ('development', 'test'):
        ds = load(root, split, scoring=True)
        row_of = {h: i for i, h in enumerate(ds.seq_ids)}
        rows = [(r.split('\t')[0], r.split('\t')[1]) for r in (ds.folder / 'annotations.tsv').read_text().splitlines()[1:]]
        first, second = ds.seq_ids[:2]
        if split == 'test':  # the first protein: curator statement and accepted review; the second: review not accepted
            kept = [r.split('\t') for r in (ds.folder / 'annotations.tsv').read_text().splitlines()]
            text = ['\t'.join(kept[0])] + ['\t'.join(r) for r in kept[1:] if r[0] != first]
            text += [f'{first}\tGO:0000002\tTAS\tGOA:X', f'{first}\tGO:0000003\treview:inference\tliterature review (this work)']
            (ds.folder / 'annotations.tsv').write_text('\n'.join(text) + '\n')
            rows = [(first, 'GO:0000002'), (first, 'GO:0000003')] + [p for p in rows if p[0] != first]
            labels, _ = ontology.label_matrix(rows, row_of)
            write('benchmark_reviews/test', [dict(protein_id=first, go_id='GO:0000003', decision='accept'),
                                             dict(protein_id=second, go_id='GO:0000002', decision='do_not_accept')])
            sp.save_npz(ds.folder / 'labels.npz', labels)
            positive = ontology.label_matrix([p for p in rows if p != (first, 'GO:0000003')], row_of)[0].toarray() > 0
            experimental = ontology.label_matrix([p for p in rows if p[0] != first], row_of)[0].toarray() > 0
        else:
            labels = ds.labels
            positive = experimental = labels.toarray() > 0
        everything = ontology.label_matrix(rows + ([(second, 'GO:0000002')] if split == 'test' else []), row_of)[0].toarray() > 0

        def table(y, mask=None):
            out = []
            for i, h in enumerate(ds.seq_ids):
                if y[i].any():
                    out.append(dict(protein_id=h, accessions=[], sequence=ds.sequences[i], cluster=ds.columns['cluster'][i],
                                    go_terms=[t for j, t in enumerate(VOCABULARY) if y[i, j]]))
                    if mask is not None:
                        out[-1]['masked_terms'] = [t for j, t in enumerate(VOCABULARY) if mask[i, j]]
            return out

        write(f'benchmark/{split}', table(labels.toarray() > 0))
        write(f'benchmark_experimental/{split}', table(experimental))
        write(f'benchmark_masked/{split}', table(positive, everything & ~positive))
        write(f'benchmark_all_goa/{split}', table(everything))
        fits[split] = root / f'{split}_scores.npy'
        np.save(fits[split], np.random.default_rng(4).random((len(ds.seq_ids), len(VOCABULARY))).astype(np.float32))
    write('benchmark_no_hit/test', [dict(protein_id=h, no_hit_up_to=10.0) for h in load(root, 'test', labels=False).seq_ids])
    return [('finch', 'archipelago', {'1': fits}), ('finch', 'hq', {'1': fits})]


@pytest.mark.parametrize('labels, proteins, pairs', [('experimental', 9, 13), ('masked', 10, 15), ('all_goa', 10, 18)])
def test_a_test_label_version_is_checked_against_the_release_and_scored(tmp_path, monkeypatch, labels, proteins, pairs):
    tiny_root(tmp_path)
    models = release_of(tmp_path)
    monkeypatch.setattr(test_labels, 'fits', lambda root, run: models)
    monkeypatch.setattr(sys, 'argv', ['test_labels', labels, '--data', str(tmp_path), '--device', 'cpu'])
    test_labels.main()
    out = tmp_path / 'results' / f'{labels}_test_labels'
    version = {r['variant']: r for r in read_table(out / 'counts.tsv')}[test_labels.VARIANT[labels]]
    assert (int(version['proteins']), int(version['propagated_pairs'])) == (proteins, pairs)
    assert {r['labels'] for r in read_table(out / 'metrics.tsv')} == {labels}


def test_a_released_label_table_that_differs_from_its_definition_stops_the_scoring(tmp_path, monkeypatch):
    tiny_root(tmp_path)
    models = release_of(tmp_path)
    path = tmp_path / 'release' / 'benchmark_all_goa' / 'test.parquet'
    rows = pq.read_table(path).to_pylist()
    rows[0]['go_terms'] = rows[0]['go_terms'][:1]
    pq.write_table(pa.Table.from_pylist(rows), path)
    monkeypatch.setattr(test_labels, 'fits', lambda root, run: models)
    monkeypatch.setattr(sys, 'argv', ['test_labels', 'all_goa', '--data', str(tmp_path), '--device', 'cpu'])
    with pytest.raises(ValueError, match='differ from their definition'):
        test_labels.main()


def test_bootstrap_weighs_proteins_by_their_cluster_draws():
    rng = np.random.default_rng(2)
    s = (rng.integers(0, 41, (30, 12)) / 40).astype(np.float32)  # ties, and scores on the thresholds
    y = rng.random((30, 12)) < 0.3
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
