"""The three simple baselines, fitted once to each training set (no seed; written as runs/<method>/<training set>/seed0/):

naive          every term scored by its frequency among the training proteins
diamond_knn    DIAMOND hits among the training proteins, bitscore-weighted share carrying each term (DeepGOPlus's rule)
embedding_knn  the nearest training proteins in ESM2-15B embedding space, similarity-weighted share carrying each term

Their raw scores are already consistent with is_a (the training labels are), so propagation leaves them unchanged.
Choices (DIAMOND-kNN's E-value cutoff, Embedding-kNN's k) are made on the development set.
"""
import numpy as np

from archipelago.evaluate import micro_auprc
from retrain.common import Ontology, fit_folder, load, write_fit

SEED = 0


class Sets:
    """A training set (with labels), the development set (with labels) and the test set (proteins only)."""

    def __init__(self, root, training_set):
        self.root, self.name = root, training_set
        self.ontology = Ontology(root)
        self.train = load(root, training_set)
        self.dev = load(root, 'development')
        self.test = load(root, 'test', labels=False)
        overlap = set(self.train.seq_ids) & (set(self.dev.seq_ids) | set(self.test.seq_ids))
        if overlap:
            raise ValueError(f'{training_set}: {len(overlap)} training proteins are development or test proteins')
        self.dev_labels = self.dev.labels.toarray()

    def dev_auprc(self, dev_scores):
        """Development micro-AUPRC of propagated scores (the scorer's headline metric). A weighted share can exceed 1 by
        a rounding error (1.0000001); it is clipped for the metric only."""
        return micro_auprc(np.clip(self.ontology.propagate(dev_scores), 0, 1), self.dev_labels)

    def write(self, method, dev_scores, test_scores, **record):
        write_fit(fit_folder(self.root, method, self.name, SEED), dev_scores, test_scores,
                  dict(method=method, training_set=self.name, seed=SEED, propagated=False,
                       dev_auprc=self.dev_auprc(dev_scores), **record))


def share(weights, labels):
    """Per query, the weighted share of its neighbours carrying each term: (weights @ labels) / row sums of weights; a query
    without neighbours scores 0. weights: sparse or dense (queries x training proteins); labels: CSR (training x terms)."""
    total = np.asarray(weights.sum(axis=1)).ravel()
    carried = weights @ labels.astype(np.float64)
    carried = carried.toarray() if hasattr(carried, 'toarray') else np.asarray(carried)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(total[:, None] > 0, carried / total[:, None], 0.0).astype(np.float32)
