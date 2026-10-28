"""Scoring against experimental-only test labels (Supplementary Information): every model of the main run scored on test
labels rebuilt from only the test annotations with an experimental evidence code.

    python -m retrain.scoring.test_labels [--data DIR] [--device cuda|cpu]

The test annotations (datasets/test/annotations.tsv) whose evidence is an experimental code (EXP, IDA, IPI, IMP, IGI, IEP,
HTP, HDA, HMP, HGI, HEP) are propagated exactly as the test labels were built (alt_id and replaced_by mapping, is_a
closure, the vocabulary); test proteins left without a label are dropped (clusters unchanged). As a check, all annotation
rows together must give the official test labels exactly. Writes results/experimental_test_labels/: the scorer's output,
contrasts.tsv (as retrain.scoring.score) and counts.tsv (proteins, propagated pairs, clusters and terms of the official
and the experimental-only labels). Needs the main run's fits.
"""
import argparse
import csv

import numpy as np

from retrain.common import Ontology, add_argument, data_root, load, log
from retrain.scoring.score import bundle, contrasts, fits, run_scorer, test_subset

EXPERIMENTAL = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--device', choices=['cuda', 'cpu'])
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    out = root / 'results' / 'experimental_test_labels'
    if out.exists():
        raise FileExistsError(f'{out} exists')
    ontology = Ontology(root)
    test = load(root, 'test', scoring=True)
    row_of = {h: i for i, h in enumerate(test.seq_ids)}
    with open(test.folder / 'annotations.tsv', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    official = test.labels.toarray() > 0
    counts = []
    for variant, kept in (('all', rows), ('exp', [r for r in rows if r['evidence'] in EXPERIMENTAL])):
        labels, dropped = ontology.label_matrix([(r['seq_id'], r['go_id']) for r in kept], row_of)
        y = labels.toarray() > 0
        if variant == 'all' and (dropped or (y != official).any()):
            raise ValueError('the test annotations do not reproduce the test labels')
        if (y & ~official).any():
            raise ValueError(f'{variant}: labels outside the official test labels')
        keep = np.flatnonzero(y.any(1))
        counts.append(dict(variant=variant, rows=len(kept), proteins=len(keep), proteins_dropped=len(row_of) - len(keep),
                           clusters=len({test.columns['cluster'][i] for i in keep}), propagated_pairs=int(y.sum()),
                           terms_with_positive=int(y.any(0).sum())))
        log(counts[-1])
    data = bundle(root)
    pairs = [(test.seq_ids[i], ontology.vocabulary[j]) for i, j in zip(*np.nonzero(y))]
    config = test_subset(root, data, root / 'scoring' / 'experimental_test_labels', {h for h, _ in pairs}, pairs,
                         fits(root, 'main'))
    run_scorer(root, config, out, a.device)
    contrasts(out)
    with open(out / 'counts.tsv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(counts[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(counts)


if __name__ == '__main__':
    main()
