"""The options for training on a published method's own split: shared proteins and vocabulary terms that are not live."""
import json
from types import SimpleNamespace

import pytest

from retrain.common import Ontology, check_disjoint

OBO = """format-version: 1.2

[Term]
id: GO:0000001
name: root
namespace: molecular_function

[Term]
id: GO:0000002
name: child
namespace: molecular_function
alt_id: GO:0000009
is_a: GO:0000001 ! root

[Term]
id: GO:0000003
name: gone
namespace: molecular_function
is_obsolete: true
"""


def ontology_root(tmp_path, vocabulary, not_live=None):
    folder = tmp_path / 'datasets' / 'ontology'
    folder.mkdir(parents=True)
    (folder / 'go-basic.obo').write_text(OBO, encoding='utf-8')
    (folder / 'vocabulary.json').write_text(json.dumps(vocabulary), encoding='utf-8')
    if not_live is not None:
        (folder / 'not_live.json').write_text(json.dumps(not_live), encoding='utf-8')
    return tmp_path


def test_terms_that_are_not_live_need_to_be_listed(tmp_path):
    vocabulary = ['GO:0000001', 'GO:0000002', 'GO:0000003', 'GO:0000009']
    with pytest.raises(ValueError):
        Ontology(ontology_root(tmp_path / 'a', vocabulary))
    with pytest.raises(ValueError):  # the list must name exactly the terms that are not live
        Ontology(ontology_root(tmp_path / 'b', vocabulary, ['GO:0000003']))
    onto = Ontology(ontology_root(tmp_path / 'c', vocabulary, ['GO:0000003', 'GO:0000009']))
    assert onto.parents == [[], [0], [], []]  # listed terms are columns without parents
    assert Ontology(ontology_root(tmp_path / 'd', vocabulary[:2])).not_live == []


def test_shared_proteins_only_when_allowed():
    training = SimpleNamespace(name='original', seq_ids=['a', 'b', 'c'])
    development, test = SimpleNamespace(seq_ids=['x']), SimpleNamespace(seq_ids=['c', 'y'])
    with pytest.raises(ValueError):
        check_disjoint(training, development, test)
    assert check_disjoint(training, development, test, allow_shared=True) == 1
    assert check_disjoint(SimpleNamespace(name='hq', seq_ids=['a']), development, test) == 0
