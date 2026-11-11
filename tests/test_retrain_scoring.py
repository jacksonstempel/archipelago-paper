"""What retrain.scoring adds to the scorer (whose end-to-end test is in the archipelago package): contrasts.tsv = the difference of the scorer's
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
from archipelago.go import information_accretion, read_parents  # noqa: E402
from archipelago.io import read_table  # noqa: E402
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


LEVELS = ('experimental', 'curated', 'computational', 'review_direct', 'review_inferred', 'review_rejected')
LEVEL_OF_EVIDENCE = {'EXP': 'experimental', 'TAS': 'curated'}


def release_of(root):
    """release/ for a tiny data root, as Hugging Face holds it (dataset v2.0), and one fit. Test proteins: the first has only
    the review's accepted assertion besides a curator statement, so it has no experimental label; the second has an
    assertion that the review did not accept."""
    (root / 'release' / 'ontology').mkdir(parents=True)
    with gzip.open(root / 'release' / 'ontology' / 'go-basic.obo.gz', 'wt') as handle:
        handle.write(OBO)

    def write(name, rows):
        path = root / 'release' / f'{name}.parquet'
        path.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.Table.from_pylist(rows), path)

    ontology = Ontology(root)
    parents = read_parents(root / 'datasets' / 'ontology' / 'go-basic.obo', VOCABULARY)
    hq = load(root, 'hq', scoring=True)
    write('go_terms/go_terms', [dict(go_id=t, column=j, parents=[VOCABULARY[p] for p in parents[j]],
                                     information_accretion=float(v))
                                for j, (t, v) in enumerate(zip(VOCABULARY, information_accretion(hq.labels.toarray() > 0, parents)))])

    def levelled(ds, rows):
        """go_terms and evidence of each protein: each propagated term with the first level whose rows give it."""
        row_of = {h: i for i, h in enumerate(ds.seq_ids)}
        level = np.zeros((len(ds.seq_ids), len(VOCABULARY)), dtype=int)
        for k, name in enumerate(LEVELS, 1):
            m = ontology.label_matrix([(h, g) for h, g, lv in rows if lv == name], row_of)[0].toarray() > 0
            level[m & (level == 0)] = k
        return [dict(go_terms=[t for j, t in enumerate(VOCABULARY) if level[i, j]],
                     evidence=[LEVELS[level[i, j] - 1] for j in range(len(VOCABULARY)) if level[i, j]])
                for i in range(len(ds.seq_ids))]

    y = hq.labels.toarray()
    write('proteins/train', [dict(protein_id=h, accessions=[], part='hq', go_terms=[t for j, t in enumerate(VOCABULARY) if y[i, j]],
                                  evidence=['experimental'] * int(y[i].sum())) for i, h in enumerate(hq.seq_ids)])
    fits = {}
    for split in ('development', 'test'):
        ds = load(root, split, scoring=True)
        row_of = {h: i for i, h in enumerate(ds.seq_ids)}
        rows = [(r.split('\t')[0], r.split('\t')[1], 'experimental')
                for r in (ds.folder / 'annotations.tsv').read_text().splitlines()[1:]]
        reviews = []
        if split == 'test':  # the first protein: curator statement and accepted review; the second: review not accepted
            first, second = ds.seq_ids[:2]
            kept = [r.split('\t') for r in (ds.folder / 'annotations.tsv').read_text().splitlines()]
            text = ['\t'.join(kept[0])] + ['\t'.join(r) for r in kept[1:] if r[0] != first]
            text += [f'{first}\tGO:0000002\tTAS\tGOA:X', f'{first}\tGO:0000003\tinference\tliterature review (this work)']
            (ds.folder / 'annotations.tsv').write_text('\n'.join(text) + '\n')
            rows = [(first, 'GO:0000002', 'curated'), (first, 'GO:0000003', 'review_inferred')] + \
                   [p for p in rows if p[0] != first] + [(second, 'GO:0000002', 'review_rejected')]
            labels, _ = ontology.label_matrix([(h, g) for h, g, lv in rows if lv != 'review_rejected'], row_of)
            sp.save_npz(ds.folder / 'labels.npz', labels)
            reviews = [dict(protein_id=first, go_id='GO:0000003', level='review_inferred', decision='accept',
                            support_type='inference'),
                       dict(protein_id=second, go_id='GO:0000002', level='review_rejected', decision='do_not_accept',
                            support_type='none_established')]
        write(f'annotations/{split}', [dict(protein_id=h, go_id=g, level=lv, decision=None, support_type=None)
                                       for h, g, lv in rows if not lv.startswith('review')] + reviews)
        write(f'proteins/{split}', [dict(protein_id=h, accessions=[], sequence=ds.sequences[i], cluster=ds.columns['cluster'][i],
                                         no_hit_up_to=10.0 if split == 'test' else None, **terms)
                                    for i, (h, terms) in enumerate(zip(ds.seq_ids, levelled(ds, rows)))])
        fits[split] = root / f'{split}_scores.npy'
        np.save(fits[split], np.random.default_rng(4).random((len(ds.seq_ids), len(VOCABULARY))).astype(np.float32))
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
    assert {r['labels'] for r in read_table(out / 'metrics.tsv') if r['split'] == 'test'} == {score.scorer_labels(labels)}

def test_the_main_test_labels_masked_are_scored_from_the_reviewed_config(tmp_path, monkeypatch):
    tiny_root(tmp_path)
    models = release_of(tmp_path)
    monkeypatch.setattr(score, 'fits', lambda root, run: models)
    for argv in (['score', 'main'], ['score', 'main', '--labels', 'masked']):
        monkeypatch.setattr(sys, 'argv', [*argv, '--data', str(tmp_path), '--device', 'cpu'])
        score.main()
    for labels, folder, pairs, left_out in (('reviewed', 'results', 16, 0), ('masked', 'results_masked', 15, 3)):
        benchmark = json.loads((tmp_path / folder / 'main' / 'inputs.json').read_text())['benchmark']
        assert (benchmark['labels'], benchmark['test_positive_pairs'], benchmark['test_masked_pairs']) == (score.scorer_labels(labels), pairs, left_out)

    def development(folder):  # the main test labels change the test set only
        return {(r['method'], r['training_set'], r['seed'], r['metric']): r['value']
                for r in read_table(tmp_path / folder / 'main' / 'metrics.tsv') if r['split'] == 'development'}

    assert development('results') == development('results_masked') and development('results')



def test_a_released_label_table_that_differs_from_its_definition_stops_the_scoring(tmp_path, monkeypatch):
    tiny_root(tmp_path)
    models = release_of(tmp_path)
    path = tmp_path / 'release' / 'proteins' / 'test.parquet'
    rows = pq.read_table(path).to_pylist()
    rows[0]['go_terms'], rows[0]['evidence'] = rows[0]['go_terms'][:1], rows[0]['evidence'][:1]
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
