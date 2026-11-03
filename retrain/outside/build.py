"""A data root for Finch on a published method's own benchmark.

    python -m retrain.outside.build {deepgose,transfew,dpfunc} --data DIR --replication REP --archipelago MAIN

REP is the method's replication root (retrain.replicate.build), MAIN the paper's data root (datasets/archipelago from the
released dataset, and Finch's embeddings of it in inputs/embeddings/<encoder>/archipelago). DIR receives:

  datasets/ontology, original_<method>, development, test   copies of REP's (the authors' GO release, vocabulary, training
                                  split, validation set, test set)
  datasets/archipelago_screened   ArchipelaGO's training proteins without a relative of the method's validation or test
                                  proteins, with ArchipelaGO's labels on the method's GO release and vocabulary:
    1. the paper's homology screen (retrain.common.screen: identical sequence, or a DIAMOND very-sensitive hit at E <= 1e-3
       at the fixed database size, either direction) against DIR's development and test sets; every protein it names goes;
    2. each asserted term mapped onto the method's GO release (alt_id -> primary id; an obsolete term with exactly one live
       replaced_by -> it; terms absent from that release dropped: Ontology.map_term), closed over is_a in that release and
       restricted to the method's vocabulary (Ontology.label_matrix); proteins left without a vocabulary label dropped.
  inputs/embeddings/<encoder>/archipelago_screened   its rows of MAIN's embeddings (Finch's encoders)
  outside/build.json      counts: screen, label mapping, proteins, terms with positives in each arm and in both
  outside/relatives.tsv   per test protein, a relative in the original training split under the same screen (identical,
                          diamond or none); ArchipelaGO has none by construction, so `none` means none in either training set
  outside/shared_terms.json   the terms with training positives in both arms

Then embed DIR's original_<method>, development and test sets with scripts/embed.py (esm2_15b_l48, prott5_l21), train
Finch on both training sets and combine (README), and score with retrain.outside.score.
"""
import argparse
import collections
import json
import os
import shutil
from pathlib import Path

import numpy as np

from retrain.common import Ontology, add_argument, data_root, load, log, read_annotations, screen, write_dataset

ARM_B = 'archipelago_screened'
ENCODERS = ('esm2_15b_l48', 'prott5_l21')


def relatives(root, proteins, work):
    """retrain.common.screen of `proteins` ({seq_id: sequence}) against the root's development and test sets: (identical
    ids, DIAMOND-hit ids, {test seq_id: 'identical' or 'diamond'} of the test proteins it links to them)."""
    exact, hits = screen(root, proteins, work, int(os.environ.get('SLURM_CPUS_PER_TASK', 8)))
    test = load(root, 'test', labels=False)
    sequences = {s.upper() for s in proteins.values()}
    found = {h: 'identical' for h, s in zip(test.seq_ids, test.sequences) if s.upper() in sequences}
    for name, column in (('source_vs_heldout.tsv', 1), ('heldout_vs_source.tsv', 0)):
        for line in open(Path(work) / name):
            held = line.split('\t')[column]
            if held.startswith('test:'):
                found.setdefault(held[5:], 'diamond')
    return exact, hits, found


def mapping_counts(ontology, pairs):
    counts = collections.Counter()
    for _, g in pairs:
        mapped = ontology.map_term(g)
        counts['dropped: absent from the release, or obsolete without one live replacement' if mapped is None else
               'same id' if mapped == g else 'alt_id -> primary id' if g in ontology.alt else 'obsolete -> replaced_by'] += 1
    return dict(counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('method', choices=['deepgose', 'transfew', 'dpfunc'])
    parser.add_argument('--replication', required=True, help="the method's replication root (retrain.replicate.build)")
    parser.add_argument('--archipelago', required=True, help="the paper's data root (datasets/archipelago and its embeddings)")
    add_argument(parser)
    a = parser.parse_args()
    root, replication, main_root = data_root(a.data), Path(a.replication), Path(a.archipelago)
    arm_a = f'original_{a.method}'
    if (root / 'datasets').exists():
        raise FileExistsError(f'{root}/datasets exists; build once')
    for name in ('ontology', arm_a, 'development', 'test'):
        shutil.copytree(replication / 'datasets' / name, root / 'datasets' / name)
    ontology = Ontology(root)
    work = root / 'outside'
    work.mkdir(parents=True)
    record = dict(method=a.method, vocabulary=len(ontology.vocabulary), not_live_terms=len(ontology.not_live))

    archipelago = load(main_root, 'archipelago', labels=False)
    proteins = dict(zip(archipelago.seq_ids, archipelago.sequences))
    exact, hits, _ = relatives(root, proteins, work / 'screen_archipelago')
    removed = exact | hits
    kept = [h for h in archipelago.seq_ids if h not in removed]
    record['archipelago'] = dict(proteins=len(proteins), identical_to_validation_or_test=len(exact), diamond_hit=len(hits),
                                 removed=len(exact | hits), kept_by_screen=len(kept))
    keep = set(kept)
    pairs = sorted({(h, g) for h, g, _, _ in read_annotations(main_root / 'datasets' / 'archipelago') if h in keep})
    record['label_mapping'] = mapping_counts(ontology, pairs)
    labels, _ = ontology.label_matrix(pairs, {h: i for i, h in enumerate(kept)})
    with_label = np.flatnonzero(labels.getnnz(axis=1) > 0)
    final, labels = [kept[i] for i in with_label], labels[with_label]
    record['archipelago'].update(without_vocabulary_label=len(kept) - len(final), training_proteins=len(final),
                                 propagated_pairs=int(labels.nnz), terms_with_positive=int((labels.getnnz(axis=0) > 0).sum()))
    kept_final = set(final)
    annotations = [(h, ontology.map_term(g), 'ArchipelaGO, mapped', f'archipelago {g}') for h, g in pairs
                   if h in kept_final and ontology.map_term(g)]
    write_dataset(root, ARM_B, final, [[] for _ in final], [proteins[h] for h in final], annotations, labels, ontology)

    original = load(root, arm_a)
    _, _, found = relatives(root, dict(zip(original.seq_ids, original.sequences)), work / 'screen_original')
    test = load(root, 'test', labels=False)
    with open(work / 'relatives.tsv', 'w') as handle:
        handle.write('seq_id\trelative_in_original\n')
        handle.writelines(f"{h}\t{found.get(h, 'none')}\n" for h in test.seq_ids)
    record['test_relatives_in_original'] = dict(collections.Counter(found.get(h, 'none') for h in test.seq_ids))
    positives_a = np.asarray(original.labels.sum(0)).ravel() > 0
    positives_b = np.asarray(labels.sum(0)).ravel() > 0
    record['terms_with_positive'] = dict(original=int(positives_a.sum()), archipelago=int(positives_b.sum()),
                                         both=int((positives_a & positives_b).sum()))
    (work / 'shared_terms.json').write_text(json.dumps([t for t, s in zip(ontology.vocabulary, positives_a & positives_b)
                                                        if s]) + '\n')

    from archipelago.io import read_matrix, write_matrix
    for encoder in ENCODERS:
        _, x = read_matrix(main_root / 'inputs' / 'embeddings' / encoder / 'archipelago', final)
        write_matrix(root / 'inputs' / 'embeddings' / encoder / ARM_B, final, np.asarray(x))
    (work / 'build.json').write_text(json.dumps(record, indent=1) + '\n')
    log(json.dumps(record))


if __name__ == '__main__':
    main()
