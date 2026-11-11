"""Scoring against the other versions of the test labels (main-text Table 2, Supplementary Information): every model of the
main run scored with the scorer's --labels.

    python -m retrain.scoring.test_labels LABELS [--data DIR] [--device cuda|cpu]

LABELS is one of
  experimental  only the test annotations with an experimental evidence code (EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP,
                HGI, HEP); test proteins left without a label are dropped (clusters unchanged)
  masked        positives: the test annotations other than the label review's; the pairs that every GOA assertion the
                review considered, accepted or not (the review rows of release/annotations/test.parquet), adds with its parent terms count
                as neither positive nor negative; every test protein is kept
  all_goa       the test annotations other than the label review's, plus every GOA assertion the review considered, accepted
                or not; every test protein is kept
-> results/<LABELS>_test_labels/

Each version is rebuilt here from the test annotations, propagated exactly as the test labels were built (alt_id and
replaced_by mapping, is_a closure, the vocabulary; a term that maps to no live term adds no label), and must equal the
released labels (the scorer's version of the same name, archipelago.dataset.labels, on release/proteins/test.parquet;
masked is the scorer's main). Further checks: all annotation rows together give the official test labels exactly, and the
review's rows among them are exactly its accepted assertions. The scorer reads the release directly; the development set
is scored on the reviewed labels, and F1 uses the threshold selected on them, as in the paper. Writes the scorer's output, contrasts.tsv (as
retrain.scoring.score) and counts.tsv (annotation rows, proteins, propagated and masked pairs, clusters and terms with a
positive, for the official labels, `all`, and the version, `exp`, `masked` or `all_goa`). Needs the main run's fits.
"""
import argparse
import csv
from collections import Counter

import numpy as np
import pyarrow.parquet as pq

from archipelago import dataset as released
from retrain.common import Ontology, add_argument, data_root, load, log, read_annotations
from retrain.scoring.score import contrasts, dataset, fits, run_scorer, scorer_labels, write_config

EXPERIMENTAL = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'}
REVIEW_SOURCE = 'literature review (this work)'  # the source of the review's accepted assertions in annotations.tsv
VARIANT = {'experimental': 'exp', 'masked': 'masked', 'all_goa': 'all_goa'}  # the name of each version in counts.tsv


def annotation_sets(root):
    """(test set, vocabulary, {name: asserted (seq_id, go_id) pairs}, propagate). The pair sets: `all` (every test
    annotation row), `exp` (the rows with an experimental evidence code), `other` (the rows other than the review's),
    `all_goa` (`other` plus every assertion the review considered, accepted or not) and `direct` (`other` plus the
    assertions the review accepted with experimental support, support_type `experimental`; the others it accepted as a
    supported inference, `inference`). propagate(pairs) gives their labels, boolean (test proteins x vocabulary). Checked:
    `all` gives the official test labels exactly, the review's rows are exactly its accepted assertions, and an assertion
    has a support type if and only if it was accepted."""
    ontology = Ontology(root)
    test = load(root, 'test', scoring=True)
    row_of = {h: i for i, h in enumerate(test.seq_ids)}

    def propagate(pairs):
        labels, dropped = ontology.label_matrix(pairs, row_of)
        if dropped:
            log(f'{len(dropped)} asserted pairs map to no live term and add no label, as in the test labels')
        return labels.toarray() > 0

    rows = read_annotations(test.folder)
    if (propagate([(h, g) for h, g, _, _ in rows]) != (test.labels.toarray() > 0)).any():
        raise ValueError('the test annotations do not reproduce the test labels')
    reviews = [r for r in pq.read_table(dataset(root) / 'annotations' / 'test.parquet',
                                        columns=['protein_id', 'go_id', 'decision', 'support_type']).to_pylist()
               if r['decision'] is not None]
    accepted = {(r['protein_id'], r['go_id']) for r in reviews if r['decision'] == 'accept'}
    if {(h, g) for h, g, _, source in rows if source == REVIEW_SOURCE} != accepted:
        raise ValueError("the review's test annotation rows are not exactly its accepted assertions")
    support = Counter((r['decision'], r['support_type']) for r in reviews)
    if any((decision == 'accept') != (kind in ('experimental', 'inference')) for decision, kind in support):
        raise ValueError(f'review decisions and support types disagree: {dict(support)}')
    log(f'label review: {len(reviews):,} assertions, {len(accepted):,} accepted; by support type {dict(support)}')
    other = [(h, g) for h, g, _, source in rows if source != REVIEW_SOURCE]
    pairs = dict(all=[(h, g) for h, g, _, _ in rows], exp=[(h, g) for h, g, e, _ in rows if e in EXPERIMENTAL],
                 other=other, all_goa=other + [(r['protein_id'], r['go_id']) for r in reviews],
                 direct=other + [(r['protein_id'], r['go_id']) for r in reviews
                                 if r['decision'] == 'accept' and r['support_type'] == 'experimental'])
    return test, ontology.vocabulary, pairs, propagate


def version_labels(pairs, propagate, labels):
    """(labels, masked pairs) of a version: boolean matrices (test proteins x vocabulary)."""
    positive = propagate(pairs[{'experimental': 'exp', 'masked': 'other', 'all_goa': 'all_goa'}[labels]])
    masked = propagate(pairs['all_goa']) & ~positive if labels == 'masked' else np.zeros_like(positive)
    return positive, masked


def write_counts(path, counts):
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(counts[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(counts)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('labels', choices=list(VARIANT))
    parser.add_argument('--device', choices=['cuda', 'cpu'])
    add_argument(parser)
    a = parser.parse_args(argv)
    root = data_root(a.data)
    out = root / 'results' / f'{a.labels}_test_labels'
    if out.exists():
        raise FileExistsError(f'{out} exists')
    test, vocabulary, pairs, propagate = annotation_sets(root)
    official = propagate(pairs['all'])
    y, masked = version_labels(pairs, propagate, a.labels)

    rows = released.proteins(dataset(root), 'test')
    if [r['protein_id'] for r in rows] != list(test.seq_ids) or released.terms(dataset(root))[0] != list(vocabulary):
        raise ValueError('the released test proteins or vocabulary differ from datasets/test and the ontology')
    y_released, masked_released = released.labels(rows, vocabulary, scorer_labels(a.labels))
    if (y_released != y).any() or (masked_released != masked).any():
        raise ValueError(f'the released {a.labels} test labels differ from their definition')

    counts = []
    for name, labels, left_out in (('all', official, np.zeros_like(official)), (VARIANT[a.labels], y, masked)):
        keep = np.flatnonzero(labels.any(1))
        counts.append(dict(variant=name, rows=len(pairs[{'masked': 'other'}.get(name, name)]), proteins=len(keep),
                           proteins_dropped=len(test.seq_ids) - len(keep),
                           clusters=len({test.columns['cluster'][i] for i in keep}), propagated_pairs=int(labels.sum()),
                           masked_pairs=int(left_out.sum()), terms_with_positive=int(labels.any(0).sum())))
        log(counts[-1])
    config = write_config(root / 'scoring' / f'{a.labels}_test_labels.json', dataset(root), fits(root, 'main'))
    run_scorer(root, config, out, a.device, labels=a.labels)
    contrasts(out)
    write_counts(out / 'counts.tsv', counts)


if __name__ == '__main__':
    main()
