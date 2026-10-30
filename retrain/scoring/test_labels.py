"""Scoring against other test labels (main-text Table 2, Supplementary Information): every model of the main run scored
on test labels rebuilt from the test annotations.

    python -m retrain.scoring.test_labels LABELS [--data DIR] [--device cuda|cpu]

LABELS is one of
  experimental  only the test annotations with an experimental evidence code (EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP,
                HGI, HEP); test proteins left without a label are dropped (clusters unchanged)
                -> results/experimental_test_labels/
  all_goa       the test annotations other than the label review's, plus every GOA assertion the review considered,
                accepted or not (release/benchmark_reviews/test.parquet); every test protein is kept
                -> results/all_goa_test_labels/

Annotations are propagated exactly as the test labels were built (alt_id and replaced_by mapping, is_a closure, the
vocabulary; a term that maps to no live term adds no label). Checks: all annotation rows together give the official test
labels exactly, and the review's rows among them are exactly its accepted assertions. Writes the scorer's output,
contrasts.tsv (as retrain.scoring.score) and counts.tsv
(annotation rows, proteins, propagated pairs, clusters and terms with a positive, for the official labels, `all`, and the
new ones, `exp` or `all_goa`). Needs the main run's fits. The labels with the reviewed assertions masked are scored by
retrain.scoring.masked.
"""
import argparse
import csv

import numpy as np
import pyarrow.parquet as pq

from retrain.common import Ontology, add_argument, data_root, load, log, read_annotations
from retrain.scoring.score import bundle, contrasts, fits, run_scorer, test_subset

EXPERIMENTAL = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'}
REVIEW_SOURCE = 'literature review (this work)'  # the source of the review's accepted assertions in annotations.tsv
FOLDER = {'exp': 'experimental_test_labels', 'all_goa': 'all_goa_test_labels'}


def annotation_sets(root):
    """(test set, vocabulary, {name: asserted (seq_id, go_id) pairs}, propagate). The pair sets: `all` (every test
    annotation row), `exp` (the rows with an experimental evidence code), `other` (the rows other than the review's) and
    `all_goa` (`other` plus every assertion the review considered, accepted or not). propagate(pairs) gives their labels,
    boolean (test proteins x vocabulary). Checked: `all` gives the official test labels exactly, and the review's rows are
    exactly its accepted assertions."""
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
    reviews = pq.read_table(root / 'release' / 'benchmark_reviews' / 'test.parquet',
                            columns=['protein_id', 'go_id', 'decision']).to_pylist()
    accepted = {(r['protein_id'], r['go_id']) for r in reviews if r['decision'] == 'accept'}
    if {(h, g) for h, g, _, source in rows if source == REVIEW_SOURCE} != accepted:
        raise ValueError("the review's test annotation rows are not exactly its accepted assertions")
    log(f'label review: {len(reviews):,} assertions, {len(accepted):,} accepted')
    other = [(h, g) for h, g, _, source in rows if source != REVIEW_SOURCE]
    pairs = dict(all=[(h, g) for h, g, _, _ in rows], exp=[(h, g) for h, g, e, _ in rows if e in EXPERIMENTAL],
                 other=other, all_goa=other + [(r['protein_id'], r['go_id']) for r in reviews])
    return test, ontology.vocabulary, pairs, propagate


def write_counts(path, counts):
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(counts[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('labels', choices=['experimental', 'all_goa'])
    parser.add_argument('--device', choices=['cuda', 'cpu'])
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    variant = {'experimental': 'exp', 'all_goa': 'all_goa'}[a.labels]
    out = root / 'results' / FOLDER[variant]
    if out.exists():
        raise FileExistsError(f'{out} exists')
    test, vocabulary, pairs, propagate = annotation_sets(root)
    official, y = propagate(pairs['all']), propagate(pairs[variant])
    if variant == 'exp' and (y & ~official).any():
        raise ValueError('experimental labels outside the official test labels')
    if variant == 'all_goa' and (official & ~y).any():
        raise ValueError('the labels with every reviewed assertion lack official test labels')
    counts = []
    for name, labels in (('all', official), (variant, y)):
        keep = np.flatnonzero(labels.any(1))
        counts.append(dict(variant=name, rows=len(pairs[name]), proteins=len(keep),
                           proteins_dropped=len(test.seq_ids) - len(keep),
                           clusters=len({test.columns['cluster'][i] for i in keep}), propagated_pairs=int(labels.sum()),
                           terms_with_positive=int(labels.any(0).sum())))
        log(counts[-1])
    data = bundle(root)
    rows = [(test.seq_ids[i], vocabulary[j]) for i, j in zip(*np.nonzero(y))]
    config = test_subset(root, data, root / 'scoring' / FOLDER[variant], {h for h, _ in rows}, rows, fits(root, 'main'))
    run_scorer(root, config, out, a.device)
    contrasts(out)
    write_counts(out / 'counts.tsv', counts)


if __name__ == '__main__':
    main()
