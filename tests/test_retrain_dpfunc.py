"""DPFunc inputs and post-processing, and the structure bank's sequence rules, on toy data (CPU, no DGL needed)."""
import hashlib
import json

import numpy as np
import pytest

from retrain.dpfunc.common import node_positions
from retrain.dpfunc.domains import scan_counts
from retrain.dpfunc.embed import window_starts
from retrain.dpfunc.inputs import edges
from retrain.structures.common import modeled_raw, read_structure, target_sequence


def test_contacts_at_10_angstrom_both_directions_in_loop_order():
    ca = np.array([[0, 0, 0], [4, 0, 0], [10, 0, 0], [20.5, 0, 0]], dtype=np.float32)
    uv = edges(ca)
    pairs = list(zip(*uv.tolist()))
    # 0-1 (4 A), 0-2 (exactly 10 A, kept), 1-2 (6 A); 2-3 is 10.5 A; no self loops; sorted by (source, target)
    assert pairs == [(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)]
    assert uv.dtype == np.int64 and edges(ca[:1]).shape == (2, 0)


def test_long_proteins_are_modeled_as_residues_1_to_1400():
    long = 'M' * 2701
    assert modeled_raw(long) == ('M' * 1400, 'residues_1_1400')
    assert modeled_raw('M' * 2700)[1] == 'full_chain'
    assert len(node_positions(2701)) == 1400 and len(node_positions(2700)) == 2700
    assert target_sequence('MUOBZXK') == 'MCKDEAK'  # proxies for residues AlphaFold2 cannot model


def test_esm_windows_cover_every_residue_with_enough_overlap():
    assert window_starts(1022) == [0]
    for length in (1023, 2000, 5000):
        starts = window_starts(length)
        assert starts[0] == 0 and starts[-1] + 1022 == length
        assert all(b - a <= 1022 - 256 for a, b in zip(starts, starts[1:]))


def test_backbone_parsing_needs_n_ca_c():
    atoms = []
    for i, name in enumerate(('MET', 'GLY'), 1):
        for atom in ('N', 'CA', 'C'):
            atoms.append(f"ATOM  {len(atoms) + 1:5d} {atom:^4s} {name} A{i:4d}    {i:8.3f}{0:8.3f}{0:8.3f}{1:6.2f}{87.5:6.2f}")
    sequence, xyz, plddt = read_structure('\n'.join(atoms).encode())
    assert sequence == 'MG' and xyz.shape == (2, 3, 3) and plddt.tolist() == [87.5, 87.5]
    with pytest.raises(ValueError):
        read_structure('\n'.join(a for a in atoms if ' C  ' not in a[11:16]).encode())


def test_interpro_counts_keep_repeated_entries_and_check_each_row(tmp_path):
    sequences = {'p1': 'MKV', 'p2': 'GG'}
    md5 = hashlib.md5(b'MKV').hexdigest()
    rows = [['p1', md5, '3', 'Pfam', 'x', 'y', '1', '3', '0.1', 'T', 'date', 'IPR000001', 'name'],
            ['p1', md5, '3', 'SMART', 'x', 'y', '1', '2', '0.1', 'T', 'date', 'IPR000001', 'name'],
            ['p1', md5, '3', 'Coils', 'x', 'y', '1', '2', '0.1', 'T', 'date']]
    path = tmp_path / 'matches.tsv'
    path.write_text('\n'.join('\t'.join(r) for r in rows) + '\n')
    assert scan_counts(sequences, path) == {'p1': {'IPR000001': 2}, 'p2': {}}
    rows[0][2] = '4'  # a match row whose sequence length differs from the protein's
    path.write_text('\n'.join('\t'.join(r) for r in rows) + '\n')
    with pytest.raises(RuntimeError):
        scan_counts(sequences, path)


OBO = """[Term]
id: GO:0000001
name: root
namespace: molecular_function

[Term]
id: GO:0000002
name: a
namespace: molecular_function
is_a: GO:0000001 ! root

[Term]
id: GO:0000003
name: b
namespace: molecular_function
relationship: part_of GO:0000002 ! a
is_a: GO:0000001 ! root
"""


def test_published_postprocessing_follows_part_of(tmp_path):
    from retrain.common import Ontology
    from retrain.dpfunc.train import published_postprocessing
    folder = tmp_path / 'datasets' / 'ontology'
    folder.mkdir(parents=True)
    (folder / 'go-basic.obo').write_text(OBO)
    (folder / 'vocabulary.json').write_text(json.dumps(['GO:0000001', 'GO:0000002', 'GO:0000003']))
    ontology = Ontology(tmp_path)
    out = published_postprocessing(ontology, np.array([[0.1, 0.2, 0.9]], dtype=np.float32))
    assert out.tolist() == [[np.float32(0.9), np.float32(0.9), np.float32(0.9)]]  # b's score reaches a through part_of
