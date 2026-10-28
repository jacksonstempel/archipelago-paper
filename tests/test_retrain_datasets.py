"""retrain.datasets on a tiny fake release: conversion, the GOA and experimental-label subsets, and the dose sets."""
import gzip
import hashlib
import json
import os
import shutil
import sys

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from retrain import common
from retrain.datasets import hq_subsets, release, selected

OBO = """format-version: 1.2

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0000001
name: a
namespace: molecular_function
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0000002
name: a child
namespace: molecular_function
is_a: GO:0000001 ! a

[Term]
id: GO:0000003
name: b
namespace: molecular_function
is_a: GO:0003674 ! molecular_function
"""
A, A2, B = 'GO:0000001', 'GO:0000002', 'GO:0000003'


def protein(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest(), sequence


P1, P2, S1, S2, S3, D1, T1 = (protein(s) for s in ('MKV', 'MKL', 'MAA', 'MCC', 'MDD', 'MEE', 'MFF'))


def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows), path)


def annotation(p, go, evidence, source, by=''):
    return dict(protein_id=p[0], go_id=go, evidence_class='', evidence=evidence, source=source, assigned_by=by)


def fake_release(root):
    r = root / 'release'
    (r / 'ontology').mkdir(parents=True)
    with gzip.open(r / 'ontology' / 'go-basic.obo.gz', 'wt') as f:
        f.write(OBO)
    write(r / 'go_terms' / 'go_terms.parquet', [dict(go_id=t, column=j) for j, t in enumerate((A, A2, B))])
    write(r / 'training' / 'train.parquet', [
        dict(protein_id=P1[0], accessions=['Q1'], sequence=P1[1], part='hq', selection_rank=None, go_terms=[A, A2]),
        dict(protein_id=P2[0], accessions=[], sequence=P2[1], part='hq', selection_rank=None, go_terms=[A, B]),
        dict(protein_id=S1[0], accessions=['X1'], sequence=S1[1], part='selected', selection_rank=1, go_terms=[B])])
    write(r / 'training_annotations' / 'train.parquet', [
        annotation(P1, A2, 'IDA', 'GOA', 'UniProt'), annotation(P2, B, 'IEA', 'GOA', 'InterPro'),
        annotation(P2, A, 'EXP', 'BRENDA'), annotation(S1, B, 'IEA', 'GOA', 'InterPro')])
    for split, p in (('development', D1), ('test', T1)):
        write(r / 'benchmark' / f'{split}.parquet',
              [dict(protein_id=p[0], accessions=[], sequence=p[1], cluster='0', go_terms=[A])])
        write(r / 'benchmark_annotations' / f'{split}.parquet', [annotation(p, A, 'IDA', 'GOA', 'UniProt')])


def run(module, monkeypatch, *args):
    monkeypatch.setattr(sys, 'argv', ['x', *map(str, args)])
    module.main()


def test_release_and_hq_subsets(tmp_path, monkeypatch):
    fake_release(tmp_path)
    run(release, monkeypatch, '--data', tmp_path)
    hq = common.load(tmp_path, 'hq')
    assert hq.seq_ids == [P1[0], P2[0]] and hq.labels.toarray().tolist() == [[1, 1, 0], [1, 0, 1]]
    arch = common.load(tmp_path, 'archipelago')
    assert arch.columns == dict(part=['hq', 'hq', 'selected'], selection_rank=['', '', '1'])
    assert common.load(tmp_path, 'test', labels=False).columns['cluster'] == ['0']

    run(hq_subsets, monkeypatch, 'goa', '--data', tmp_path)
    goa = common.load(tmp_path, 'goa')  # P2 has only a computational GOA annotation: not a GOA member
    assert goa.seq_ids == [P1[0]] and goa.labels.toarray().tolist() == [[1, 1, 0]]
    run(hq_subsets, monkeypatch, 'hq_explabels', '--data', tmp_path)
    assert common.load(tmp_path, 'hq_explabels').labels.toarray().tolist() == [[1, 1, 0], [1, 0, 0]]


def test_hq_subsets_follow_the_support_file(tmp_path, monkeypatch):
    """With release/hq_support.tsv, only admitted GOA annotations make a protein a GOA member, and rule completions come
    from the file; the paper's sizes are enforced."""
    fake_release(tmp_path)
    run(release, monkeypatch, '--data', tmp_path)
    (tmp_path / 'release' / 'hq_support.tsv').write_text(
        'seq_id\tgo_id\tevidence\tassigned_by\troute\n'
        f'{P2[0]}\t{B}\tIDA\tUniProt\tadmitted_goa\n'  # P2 admitted, P1's IDA row not
        f'{P1[0]}\t{A2}\t\t\trule_completion\n', encoding='utf-8')
    monkeypatch.setattr(hq_subsets, 'PAPER', {'goa': (1, 1), 'hq_explabels': (2, 4)})
    run(hq_subsets, monkeypatch, 'goa', '--data', tmp_path)
    goa = common.load(tmp_path, 'goa')
    assert goa.seq_ids == [P2[0]] and goa.labels.toarray().tolist() == [[0, 0, 1]]
    run(hq_subsets, monkeypatch, 'hq_explabels', '--data', tmp_path)
    assert common.load(tmp_path, 'hq_explabels').labels.toarray().tolist() == [[1, 1, 0], [1, 0, 1]]
    shutil.rmtree(tmp_path / 'datasets' / 'goa')
    monkeypatch.setattr(hq_subsets, 'PAPER', {'goa': (1, 2), 'hq_explabels': (2, 4)})
    with pytest.raises(ValueError, match="differs from the paper's"):
        run(hq_subsets, monkeypatch, 'goa', '--data', tmp_path)


def test_candidates_from_release(tmp_path, monkeypatch):
    from retrain.datasets import candidates
    fake_release(tmp_path)
    run(release, monkeypatch, '--data', tmp_path)
    write(tmp_path / 'release' / 'candidates' / 'proteins.parquet', [
        dict(protein_id=S2[0], accessions=['C2'], sequence=S2[1], go_terms=[A, A2]),
        dict(protein_id=S3[0], accessions=[], sequence=S3[1], go_terms=[B])])
    rows = [dict(protein_id=S2[0], go_id=A2, evidence='IEA', assigned_by='InterPro'),
            dict(protein_id=S3[0], go_id=B, evidence='IEA', assigned_by='InterPro')]
    write(tmp_path / 'release' / 'candidates' / 'annotations.parquet', rows)
    run(candidates, monkeypatch, '--data', tmp_path)
    ds = common.load(tmp_path, 'candidates')
    assert ds.seq_ids == [S2[0], S3[0]] and ds.labels.toarray().tolist() == [[1, 1, 0], [0, 0, 1]]
    shutil.rmtree(tmp_path / 'datasets' / 'candidates')
    write(tmp_path / 'release' / 'candidates' / 'annotations.parquet', rows[:1])  # S3's label now unexplained
    with pytest.raises(ValueError, match='differ from the released'):
        run(candidates, monkeypatch, '--data', tmp_path)


def test_hq_random_from_the_released_draw(tmp_path, monkeypatch):
    from retrain.datasets import random_draw
    fake_release(tmp_path)
    run(release, monkeypatch, '--data', tmp_path)
    write(tmp_path / 'release' / 'hq_random' / 'proteins.parquet', [
        dict(protein_id=S3[0], accessions=['R3'], sequence=S3[1], go_terms=[B]),
        dict(protein_id=S2[0], accessions=[], sequence=S2[1], go_terms=[A, A2])])
    write(tmp_path / 'release' / 'hq_random' / 'annotations.parquet', [
        dict(protein_id=S3[0], go_id=B, evidence='IEA', assigned_by='InterPro'),
        dict(protein_id=S2[0], go_id=A2, evidence='IEA', assigned_by='')])
    monkeypatch.setattr(random_draw, 'KEEP', 2)
    run(random_draw, monkeypatch, '--data', tmp_path)
    ds = common.load(tmp_path, 'hq_random')
    assert ds.seq_ids == [P1[0], P2[0], S3[0], S2[0]] and ds.columns['order_rank'] == ['', '', '1', '2']
    assert ds.labels.toarray()[2:].tolist() == [[0, 0, 1], [1, 1, 0]]
    assert (S2[0], A2, 'IEA', 'GOA') in common.read_annotations(ds.folder)


def test_release_refuses_labels_its_annotations_do_not_give(tmp_path, monkeypatch):
    fake_release(tmp_path)
    table = pq.read_table(tmp_path / 'release' / 'training' / 'train.parquet').to_pylist()
    table[1]['go_terms'] = [A, A2, B]  # A2 is not implied by P2's annotations
    write(tmp_path / 'release' / 'training' / 'train.parquet', table)
    with pytest.raises(ValueError, match='differ from the released'):
        run(release, monkeypatch, '--data', tmp_path)


def test_dose_sets_extend_archipelago(tmp_path, monkeypatch):
    fake_release(tmp_path)
    run(release, monkeypatch, '--data', tmp_path)
    ontology = common.Ontology(tmp_path)
    ids, seqs = [S1[0], S2[0], S3[0]], [S1[1], S2[1], S3[1]]
    labels, _ = ontology.label_matrix([(S1[0], B), (S2[0], A2), (S3[0], B)], {h: i for i, h in enumerate(ids)})
    common.write_dataset(tmp_path, 'candidates', ids, [[], [], []], seqs,
                         [(S1[0], B, 'IEA', 'GOA:InterPro'), (S2[0], A2, 'IEA', 'GOA:InterPro'), (S3[0], B, 'IEA', 'GOA:X')],
                         labels, ontology)
    (tmp_path / 'selection').mkdir()
    (tmp_path / 'selection' / 'coverage_ranking.tsv').write_text(
        'rank\tseq_id\tscore\n' + ''.join(f'{k}\t{h}\t0.5\n' for k, h in enumerate(ids, 1)))
    monkeypatch.setattr(selected, 'DOSES', {'archipelago_3': 3})
    run(selected, monkeypatch, 'dose', '--data', tmp_path)
    dose = common.load(tmp_path, 'archipelago_3')
    assert dose.seq_ids == [P1[0], P2[0], *ids]
    assert dose.labels.toarray()[3].tolist() == [1, 1, 0]
    assert json.loads((tmp_path / 'datasets' / 'ontology' / 'vocabulary.json').read_text(encoding='utf-8')) == [A, A2, B]
    assert np.array_equal(dose.labels.toarray()[:3], common.load(tmp_path, 'archipelago').labels.toarray())


@pytest.mark.skipif(not (os.environ.get('DIAMOND') or shutil.which('diamond')), reason='needs DIAMOND for the screen')
def test_control_sets_add_screened_candidates_only(tmp_path, monkeypatch):
    """Each Fig. 4b control set is HQ plus exactly KEEP candidates, none of them an HQ or benchmark protein, and a candidate
    identical to a test protein is screened out."""
    fake_release(tmp_path)
    run(release, monkeypatch, '--data', tmp_path)
    ontology = common.Ontology(tmp_path)
    extra = [protein(s) for s in ('MGG', 'MHH')]
    candidates = [S1, S2, T1, S3, *extra]  # T1: the sequence of the test protein, planted among the candidates
    ids = [p[0] for p in candidates]
    labels, _ = ontology.label_matrix([(h, B) for h in ids], {h: i for i, h in enumerate(ids)})
    common.write_dataset(tmp_path, 'candidates', ids, [[] for _ in ids], [p[1] for p in candidates],
                         [(h, B, 'IEA', 'GOA:InterPro') for h in ids], labels, ontology)
    (tmp_path / 'selection').mkdir()
    (tmp_path / 'selection' / 'coverage_ranking_unfiltered.tsv').write_text(
        'rank\tseq_id\tscore\n' + ''.join(f'{k}\t{h}\t0.5\n' for k, h in enumerate(ids, 1)))
    monkeypatch.setattr(selected, 'KEEP', 3)
    monkeypatch.setattr(selected, 'DRAW', len(ids))
    monkeypatch.setattr(selected, 'POOL', len(ids))
    monkeypatch.setattr(selected, 'ELIGIBLE', len(ids))  # B is held by one HQ protein: every candidate is rare-function
    run(selected, monkeypatch, 'controls', '--data', tmp_path, '--threads', 1)
    hq = set(common.load(tmp_path, 'hq', labels=False).seq_ids)
    for name in ('pool_random', 'filter_random', 'coverage_unfiltered'):
        ds = common.load(tmp_path, name)
        added = ds.seq_ids[len(hq):]
        assert set(ds.seq_ids[:len(hq)]) == hq and len(added) == 3
        assert T1[0] not in added and not set(added) & (hq | {D1[0]})
