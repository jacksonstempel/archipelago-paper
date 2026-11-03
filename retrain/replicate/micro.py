"""The paper's micro-AUPRC on a published method's own benchmark.

archipelago.evaluate.micro_auprc (average precision over all protein-term pairs pooled, over every distinct score) on the
propagated prediction and truth matrices that the method's own evaluation code scores, restricted to its vocabulary
without the MF root:
  DeepGO-SE  predictions propagated over go.obo with every relation (the models' scores as its pipeline propagates them),
             truth = the records' prop_annotations
  TransFew   CAFA-evaluator's matrices: predict.py's terms with scores > 0.01 at full precision (a prediction file as written,
             the authors', at its three decimals), max-propagated over TFewData's go-basic.obo (is_a, part_of), and its
             ground truth propagated the same way
  DPFunc     new_compute_performance_deepgoplus's max-propagated predictions and closed truth, as it hands them to its fmax
Proteins without a positive among the evaluated terms are left out. The papers' own protein-centric trapezoid AUPR is kept
by the scorers as an option (--paper-aupr): the check that the papers' protocols are reproduced.
"""
import numpy as np

ROOT = 'GO:0003674'


def micro(scores, truth):
    from archipelago.evaluate import micro_auprc
    truth = np.asarray(truth, dtype=bool)
    keep = truth.any(1)
    return micro_auprc(np.asarray(scores, dtype=np.float32)[keep], truth[keep])


def deepgose(pred, prop_annotations, terms, go_set):
    """pred: propagated scores (records x terms); prop_annotations: each record's truth; go_set: the MF terms of go.obo
    without the root."""
    cols = [j for j, t in enumerate(terms) if t in go_set]
    index = {terms[j]: k for k, j in enumerate(cols)}
    truth = np.zeros((len(pred), len(cols)), dtype=bool)
    for i, annotations in enumerate(prop_annotations):
        for t in annotations:
            k = index.get(t)
            if k is not None:
                truth[i, k] = True
    return micro(np.asarray(pred)[:, cols], truth)


class TransFew:
    """CAFA-evaluator's matrices of TransFew's evaluation (TFewData/obo/go-basic.obo, is_a and part_of, max)."""

    def __init__(self, tfew, vocabulary):
        import cafaeval.parser as cp
        self.cp = cp
        self.ontology = cp.obo_parser(str(tfew / 'obo' / 'go-basic.obo'), ('is_a', 'part_of'),
                                      str(tfew / 'test' / 'output_t1_t2' / 'ia_prop.txt'), True)
        self.mf = self.ontology['molecular_function']
        self.position = {t: d['index'] for t, d in self.mf.terms_dict.items()}
        self.vocabulary = [t for t in vocabulary if t != ROOT and t in self.position]

    def score(self, gt_file, keep=None, scores=None, proteins=None, terms=None, columns=None, pred_file=None):
        """micro-AUPRC against gt_file: predictions from `scores` ({seq_id: row}; predict.py's `terms` at `columns`, > 0.01)
        of `proteins` [(accession, seq_id)], or from a prediction file as written; columns restricted to `keep`."""
        gt = self.cp.gt_parser(str(gt_file), self.ontology)['molecular_function']
        if pred_file is not None:
            m = self.cp.pred_parser(str(pred_file), self.ontology, {'molecular_function': gt}, 'max')['molecular_function'].matrix
        else:
            m = np.zeros(gt.matrix.shape)
            target, columns, seq = np.array([self.position[t] for t in terms]), np.asarray(columns), dict(proteins)
            for accession, i in gt.ids.items():
                values = np.asarray(scores[seq[accession]], dtype=float)[columns]
                hit = values > 0.01
                m[i, target[hit]] = values[hit]
            self.cp.propagate(m, self.mf, self.mf.order, mode='max')
        cols = [self.position[t] for t in self.vocabulary if keep is None or t in keep]
        return micro(m[:, cols], gt.matrix[:, cols])


def dpfunc(evaluation, frame, obo, vocabulary, keeps=(None,)):
    """new_compute_performance_deepgoplus on the frame, keeping the matrices it hands to its fmax: (its result (fmax, aupr,
    threshold), [micro-AUPRC over the vocabulary without the root, restricted to each of `keeps`])."""
    captured, fmax, terms_of = {}, evaluation.fmax, evaluation.Ontology.get_namespace_terms

    def take_fmax(targets, scores):
        captured['truth'], captured['scores'] = np.asarray(targets) > 0, np.asarray(scores)
        return fmax(targets, scores)

    def take_terms(self, namespace):
        captured['terms'] = terms_of(self, namespace)  # the set the function lists its columns from
        return captured['terms']
    evaluation.fmax, evaluation.Ontology.get_namespace_terms = take_fmax, take_terms
    try:
        result = evaluation.new_compute_performance_deepgoplus(frame, obo, 'mf')
    finally:
        evaluation.fmax, evaluation.Ontology.get_namespace_terms = fmax, terms_of
    where = {t: j for j, t in enumerate(captured['terms'])}
    out = []
    for keep in keeps:
        cols = [where[t] for t in vocabulary if t != ROOT and t in where and (keep is None or t in keep)]
        out.append(micro(np.clip(captured['scores'][:, cols], 0, 1), captured['truth'][:, cols]))
    return result, out
