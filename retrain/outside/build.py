"""A data root for Finch on a published method's own benchmark.

    python -m retrain.outside.build {deepgose,transfew,dpfunc} --data DIR --replication REP --archipelago MAIN [--addback NAME]
        [--screen-original]

REP is the method's replication root (retrain.replicate.build), MAIN the paper's data root (datasets/archipelago from the
released dataset, and Finch's embeddings of it in inputs/embeddings/<encoder>/archipelago). DIR receives:

  datasets/ontology, original_<method>, development, test   copies of REP's (the authors' GO release, vocabulary, training
                                  split, validation set, test set)
  datasets/archipelago_screened   ArchipelaGO's training proteins (and, with --addback, the proteins of MAIN's dataset
                                  NAME) without a relative of the method's test proteins, with their labels on the method's
                                  GO release and vocabulary:
    1. the paper's homology screen (retrain.common.screen: identical sequence, or a DIAMOND very-sensitive hit at E <= 1e-3
       at the fixed database size, either direction) against DIR's test set (the validation set only selects Finch's
       checkpoint and F1 threshold, as for the original split); every protein it names goes;
    2. each asserted term mapped onto the method's GO release (alt_id -> primary id; an obsolete term with exactly one live
       replaced_by -> it; terms absent from that release dropped: Ontology.map_term), closed over the relations of the
       method's own training labels in that release (RELATIONS: DeepGO-SE is_a and every relationship, as deepgo2's
       Ontology with_rels; TransFew and DPFunc is_a and part_of; each method's training labels are closed under exactly
       these) and restricted to the method's vocabulary; proteins left without a vocabulary label dropped. Column `part`
       says where each protein comes from.
  inputs/embeddings/<encoder>/archipelago_screened   its rows of MAIN's embeddings (Finch's encoders)
  outside/build.json      counts: screen, label mapping, proteins (per part), terms with positives in each arm and in both
  outside/relatives.tsv   per test protein, a relative in the original training split under the same screen (identical,
                          diamond or none); arm B has none by construction, so `none` means none in either training set
  outside/shared_terms.json   the terms with training positives in both arms

--screen-original: also datasets/original_<method>_screened, the authors' training split without the proteins the same
screen links to the test set (counts in build.json); train Finch on it as on original_<method>, score it with
retrain.outside.score --screened-original (the paper: TransFew and DPFunc).

--addback NAME: the paper's runs added back the proteins ArchipelaGO left out only because of the paper's own benchmark
(its 2,500 development and 5,000 test proteins and the 4,063 HQ proteins related to them), labelled as HQ training
proteins are labelled (not with the reviewed benchmark labels); MAIN must hold them as datasets/NAME with their
embeddings (scripts/embed.py --set NAME).

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
import scipy.sparse as sp

from retrain.common import Ontology, add_argument, data_root, load, log, read_annotations, screen, write_dataset

ARM_B = 'archipelago_screened'
ARM_A = 'original_{}'
ARM_A_SCREENED = 'original_{}_screened'  # --screen-original: arm A without the test proteins' relatives
ENCODERS = ('esm2_15b_l48', 'prott5_l21')
RELATIONS = {'deepgose': 'all', 'transfew': ('part_of',), 'dpfunc': ('part_of',)}  # besides is_a


def relatives(root, proteins, work):
    """retrain.common.screen of `proteins` ({seq_id: sequence}) against the root's test set: (identical ids, DIAMOND-hit
    ids, {test seq_id: 'identical' or 'diamond'} of the test proteins it links to them)."""
    exact, hits = screen(root, proteins, work, int(os.environ.get('SLURM_CPUS_PER_TASK', 8)), against=('test',))
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


def ancestors_over(obo, relations):
    """A function term -> {term and its ancestors} over is_a and `relations` ('all': every relationship type) in the whole
    OBO file (all namespaces); breadth-first, so cycles through has_part / part_of are harmless."""
    parents, term = collections.defaultdict(list), None
    for line in open(obo, encoding='utf-8-sig'):
        line = line.rstrip('\n')
        if line.startswith('['):
            term = '' if line == '[Term]' else None
        elif term is not None and line.startswith('id: '):
            term = line[4:].strip()
        elif term and line.startswith('is_a: '):
            parents[term].append(line[6:].split()[0])
        elif term and line.startswith('relationship: '):
            kind, target = line[len('relationship: '):].split()[:2]
            if relations == 'all' or kind in relations:
                parents[term].append(target)
    memo = {}

    def ancestors(t):
        if t not in memo:
            seen, frontier = {t}, [t]
            while frontier:
                frontier = [p for x in frontier for p in parents.get(x, ()) if p not in seen]
                seen.update(frontier)
            memo[t] = frozenset(seen)
        return memo[t]
    return ancestors


def label_matrix(ontology, ancestors, pairs, row_of):
    """Asserted (key, go_id) pairs -> CSR uint8 labels: mapped (Ontology.map_term), closed over `ancestors`, restricted to
    the vocabulary."""
    rows, cols = [], []
    for key, g in pairs:
        t = ontology.map_term(g)
        if t is not None:
            for x in ancestors(t):
                j = ontology.index.get(x)
                if j is not None:
                    rows.append(row_of[key])
                    cols.append(j)
    m = sp.coo_matrix((np.ones(len(rows), np.int32), (rows, cols)), shape=(len(row_of), len(ontology.vocabulary))).tocsr()
    m.data[:] = 1
    return m.astype(np.uint8)


def screened_original(root, method, arm, exact, hits):
    """Arm A screened: the authors' training proteins without a relative of the method's test proteins under arm B's
    screen (identical sequence, or a DIAMOND very-sensitive hit, test set only); labels and annotations as arm A. Writes
    datasets/original_<method>_screened; returns its counts."""
    from retrain.common import Ontology, read_annotations, write_dataset
    removed = (exact | hits) & set(arm.seq_ids)
    keep = [i for i, h in enumerate(arm.seq_ids) if h not in removed]
    ids = [arm.seq_ids[i] for i in keep]
    kept = set(ids)
    annotations = [r for r in read_annotations(root / 'datasets' / ARM_A.format(method)) if r[0] in kept]
    write_dataset(root, ARM_A_SCREENED.format(method), ids, [arm.accessions[i].split(',') if arm.accessions[i] else []
                                                             for i in keep],
                  [arm.sequences[i] for i in keep], annotations, arm.labels[keep], Ontology(root),
                  columns={k: [v[i] for i in keep] for k, v in arm.columns.items()})
    return dict(training_proteins=len(arm.seq_ids), removed=len(removed), identical_to_test=len(exact & removed),
                diamond_hit_only=len(removed - exact), kept=len(ids))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('method', choices=['deepgose', 'transfew', 'dpfunc'])
    parser.add_argument('--replication', required=True, help="the method's replication root (retrain.replicate.build)")
    parser.add_argument('--archipelago', required=True, help="the paper's data root (datasets/archipelago and its embeddings)")
    parser.add_argument('--addback', help="a dataset of MAIN added to ArchipelaGO's training proteins before the screen")
    parser.add_argument('--screen-original', action='store_true',
                        help='also datasets/original_<method>_screened: the training split without the relatives of the test '
                             'proteins under the same screen as archipelago_screened (the paper: TransFew and DPFunc)')
    add_argument(parser)
    a = parser.parse_args()
    root, replication, main_root = data_root(a.data), Path(a.replication), Path(a.archipelago)
    arm_a = f'original_{a.method}'
    if (root / 'datasets').exists():
        raise FileExistsError(f'{root}/datasets exists; build once')
    for name in ('ontology', arm_a, 'development', 'test'):
        shutil.copytree(replication / 'datasets' / name, root / 'datasets' / name)
    ontology = Ontology(root)
    ancestors = ancestors_over(root / 'datasets' / 'ontology' / 'go-basic.obo', RELATIONS[a.method])
    work = root / 'outside'
    work.mkdir(parents=True)
    record = dict(method=a.method, vocabulary=len(ontology.vocabulary), not_live_terms=len(ontology.not_live),
                  relations=['is_a', RELATIONS[a.method]])

    proteins, part, pairs = {}, {}, set()
    for name in ('archipelago', *([a.addback] if a.addback else [])):
        ds = load(main_root, name, labels=False)
        if set(ds.seq_ids) & set(proteins):
            raise ValueError(f'{name} shares proteins with ArchipelaGO')
        proteins.update(zip(ds.seq_ids, ds.sequences))
        part.update(dict.fromkeys(ds.seq_ids, name))
        pairs |= {(h, g) for h, g, _, _ in read_annotations(main_root / 'datasets' / name)}
    candidates = list(proteins)
    exact, hits, _ = relatives(root, proteins, work / 'screen_archipelago')
    removed = exact | hits
    kept = [h for h in candidates if h not in removed]
    keep = set(kept)
    pairs = sorted((h, g) for h, g in pairs if h in keep)
    labels = label_matrix(ontology, ancestors, pairs, {h: i for i, h in enumerate(kept)})
    with_label = np.flatnonzero(labels.getnnz(axis=1) > 0)
    final, labels = [kept[i] for i in with_label], labels[with_label]
    count = lambda hs: dict(collections.Counter(part[h] for h in hs))  # noqa: E731
    record['arm_b'] = dict(candidates=count(candidates), identical_to_test=count(exact), removed_by_screen=count(removed),
                           kept_by_screen=count(kept), without_vocabulary_label=count(keep - set(final)),
                           proteins=count(final), total=len(final), propagated_pairs=int(labels.nnz),
                           terms_with_positive=int((labels.getnnz(axis=0) > 0).sum()))
    record['label_mapping'] = {name: mapping_counts(ontology, [(h, g) for h, g in pairs if part[h] == name])
                               for name in sorted(set(part.values()))}
    kept_final = set(final)
    annotations = [(h, ontology.map_term(g), 'asserted, mapped onto this GO release', f'{part[h]} {g}') for h, g in pairs
                   if h in kept_final and ontology.map_term(g)]
    write_dataset(root, ARM_B, final, [[] for _ in final], [proteins[h] for h in final], annotations, labels, ontology,
                  columns={'part': [part[h] for h in final]})

    original = load(root, arm_a)
    exact_a, hits_a, found = relatives(root, dict(zip(original.seq_ids, original.sequences)), work / 'screen_original')
    if a.screen_original:
        record['original_screened'] = screened_original(root, a.method, original, exact_a, hits_a)
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
        rows = {}
        for name in sorted(set(part[h] for h in final)):
            ids = [h for h in final if part[h] == name]
            _, x = read_matrix(main_root / 'inputs' / 'embeddings' / encoder / name, ids)
            rows.update(zip(ids, np.asarray(x)))
        write_matrix(root / 'inputs' / 'embeddings' / encoder / ARM_B, final, np.stack([rows[h] for h in final]))
    (work / 'build.json').write_text(json.dumps(record, indent=1) + '\n')
    log(json.dumps(record))


if __name__ == '__main__':
    main()
