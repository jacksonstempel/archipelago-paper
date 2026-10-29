"""scripts/export_package.py on CPU: a tiny Finch package (tiny readouts, 5 GO terms, fake language models) gets the
development micro-Fmax threshold of its own seed ensemble, and checkpoints of another vocabulary are refused."""
import json
import zlib
from pathlib import Path

import numpy as np
import pytest
import torch

from archipelago import Finch, embed, finch, readout
from archipelago.evaluate import micro_fmax
from archipelago.go import Propagation, read_parents
from archipelago.io import write_matrix
from scripts import export_package

OBO = """format-version: 1.2
data-version: releases/test

[Term]
id: GO:0000001
name: root

[Term]
id: GO:0000002
name: a
is_a: GO:0000001 ! root

[Term]
id: GO:0000003
name: b
is_a: GO:0000001 ! root

[Term]
id: GO:0000004
name: child of a and b
is_a: GO:0000002 ! a
is_a: GO:0000003 ! b

[Term]
id: GO:0000005
name: grandchild
is_a: GO:0000004 ! child of a and b
is_a: GO:0000099 ! outside the vocabulary
"""
VOCABULARY = ["GO:0000001", "GO:0000002", "GO:0000003", "GO:0000004", "GO:0000005"]
ANCESTORS = [[0], [0, 1], [0, 2], [0, 1, 2, 3], [0, 1, 2, 3, 4]]
DIMS = {"esm2_15b_l48": 6, "prott5_l21": 4}
SEEDS = [11, 12]
WIDTH = 8
SEQ_A, SEQ_B = "MKTAYIAKQRQISFVKSHFSRQ", "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQ"


def vector(name, sequence):
    """A fake language model's embedding: fixed per model and sequence."""
    return np.random.default_rng(zlib.crc32(f"{name}:{sequence}".encode())).normal(size=DIMS[name]).astype(np.float32)


class FakeLanguageModel:
    alive, calls = 0, []

    def __init__(self, name):
        assert FakeLanguageModel.alive == 0, "two language models in memory at once"
        FakeLanguageModel.alive += 1
        self.name = name

    def __call__(self, sequences):
        FakeLanguageModel.calls.append((self.name, list(sequences)))
        return np.stack([vector(self.name, s) for s in sequences])

    def __del__(self):
        FakeLanguageModel.alive -= 1


@pytest.fixture(autouse=True)
def fake_language_models(monkeypatch):
    for name in DIMS:
        monkeypatch.setitem(embed.LOADERS, name, lambda device, name=name: FakeLanguageModel(name))
    FakeLanguageModel.calls = []


def write_ontology(root, vocabulary):
    ontology = root / "datasets/ontology"
    ontology.mkdir(parents=True)
    (ontology / "go-basic.obo").write_text(OBO, encoding="utf-8")
    (ontology / "vocabulary.json").write_text(json.dumps(vocabulary), encoding="utf-8")


@pytest.fixture(scope="module")
def package(tmp_path_factory):
    import scipy.sparse as sp
    root, rng = tmp_path_factory.mktemp("finch"), np.random.default_rng(0)
    write_ontology(root, VOCABULARY)
    config = json.loads((Path(__file__).resolve().parents[1] / "configs/finch.json").read_text(encoding="utf-8"))
    config["seeds"] = config["release_seeds"] = SEEDS
    for name in DIMS:
        config["readouts"][name]["width"] = WIDTH
    (root / "config.json").write_text(json.dumps(config), encoding="utf-8")

    for name, dim in DIMS.items():  # checkpoints in retrain.finch's format
        for seed in SEEDS:
            torch.manual_seed(seed + dim)
            model = readout.Readout(dim, WIDTH, len(VOCABULARY), 0.7)
            with torch.no_grad():
                model.network[4].weight.mul_(8.0)  # spread the scores
            folder = root / f"runs/finch/archipelago/seed{seed}/{name}"
            folder.mkdir(parents=True)
            torch.save(dict(state_dict=model.state_dict(), readout=name, seed=seed, recipe=config["readouts"][name],
                            input_dim=dim, terms=len(VOCABULARY), vocabulary=VOCABULARY, best_step=1000),
                       folder / "readout.pt")

    ids = [f"dev{i}" for i in range(30)]
    sequences = ["".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), 20)) for _ in ids]
    labels = np.zeros((len(ids), len(VOCABULARY)), dtype=np.uint8)
    for i in range(len(ids)):
        labels[i, ANCESTORS[rng.integers(0, len(VOCABULARY))]] = 1
    development = root / "datasets/development"
    development.mkdir()
    (development / "proteins.tsv").write_text("seq_id\taccessions\n" + "".join(f"{p}\t\n" for p in ids), encoding="utf-8")
    (development / "sequences.fasta").write_text("".join(f">{p}\n{s}\n" for p, s in zip(ids, sequences)),
                                                 encoding="utf-8")
    sp.save_npz(development / "labels.npz", sp.csr_matrix(labels))
    for name in DIMS:
        write_matrix(root / "inputs/embeddings" / name / "development", ids,
                     np.stack([vector(name, s) for s in sequences]).astype(np.float16))

    export_package.main(["--data", str(root), "--config", str(root / "config.json"), "--device", "cpu",
                         "--out", str(root / "package")])
    return dict(root=root, config=config, dev_sequences=sequences, dev_labels=labels)


def expected_scores(package, sequences):
    """Finch's seed ensemble computed directly from the training checkpoints, as the scorer's ensemble: per seed,
    readout.predict, finch.combine and go.Propagation; then the mean of the seeds' propagated scores."""
    propagate = Propagation(read_parents(package["root"] / "datasets/ontology/go-basic.obo", VOCABULARY))
    per_seed = []
    for seed in SEEDS:
        raw = {}
        for name in DIMS:
            saved = torch.load(package["root"] / f"runs/finch/archipelago/seed{seed}/{name}/readout.pt",
                               weights_only=True)
            model = readout.Readout(DIMS[name], WIDTH, len(VOCABULARY), 0.7)
            model.load_state_dict(saved["state_dict"])
            raw[name] = readout.predict(model, np.stack([vector(name, s) for s in sequences]).astype(np.float16))
        per_seed.append(propagate(finch.combine(raw, package["config"]["weights"])))
    return np.mean(np.array(per_seed, dtype=np.float64), 0).astype(np.float32)


def test_package_threshold_is_the_development_micro_fmax_threshold(package):
    model = Finch.from_pretrained(str(package["root"] / "package"), device="cpu")
    assert model.threshold == micro_fmax(expected_scores(package, package["dev_sequences"]), package["dev_labels"])[1]
    assert model.names == ["root", "a", "b", "child of a and b", "grandchild"]


def test_export_rejects_checkpoints_of_another_vocabulary(package, tmp_path):
    write_ontology(tmp_path, [VOCABULARY[1], VOCABULARY[0]] + VOCABULARY[2:])  # same width, other column meanings
    with pytest.raises(ValueError, match="vocabulary"):
        export_package.main(["--data", str(tmp_path), "--config", str(package["root"] / "config.json"),
                             "--device", "cpu", "--out", str(tmp_path / "package"),
                             "--checkpoints", str(package["root"]) + "/runs/finch/archipelago/seed{seed}/{readout}/readout.pt"])
