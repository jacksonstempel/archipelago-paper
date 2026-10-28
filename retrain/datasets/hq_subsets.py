"""Two training sets made from the HQ proteins and their annotations (datasets/hq):

    python -m retrain.datasets.hq_subsets goa [--data DIR]
    python -m retrain.datasets.hq_subsets hq_explabels [--data DIR]

goa           The HQ proteins with GOA support: a GOA annotation with an experimental evidence code (EXP, IDA, IPI, IMP,
              IGI, IEP, HTP, HDA, HMP, HGI, HEP) or a curator statement (TAS, IC) that the HQ build admitted. Labels: all
              their GOA annotations, any evidence code (including the older IBA/IEA panel labels); none of the
              specialist-database, UniProt catalytic-activity, label-completion-rule or review labels. Rows in HQ order;
              column `experimental` (1 if an admitted GOA annotation has an experimental code, 0 if only TAS/IC). The
              first rung of the paper's ladder: 69,412 proteins, 911,848 labels.
hq_explabels  The same 87,797 proteins as HQ, labelled only with what makes a protein HQ: the admitted GOA experimental and
              TAS/IC annotations, specialist databases, UniProt catalytic activities and the label-completion rules; every
              other (computational) GOA annotation is dropped. Every label is an hq label (checked). 644,115 labels.

Exact sets need one release file, release/hq_support.tsv (to be released with the dataset; about 150,000 rows, 12 MB),
tab-separated with a header:
  seq_id       an HQ protein
  go_id        the GO term as asserted
  evidence     the evidence code (empty for rule completions)
  assigned_by  the GOA assigning group (empty for rule completions)
  route        admitted_goa: a GOA experimental or TAS/IC annotation that the HQ build admitted as experimental support
               (isoform-specific accession, object type protein, RBP2GO excluded, conflict list applied);
               rule_completion: a label added by the label-completion rules (propagated pairs)
With it, both sets are rebuilt exactly and their sizes must be the paper's (checked).

Without it, the build falls back, with a warning, to the released annotations alone, which differ in two ways: 200 GOA
pairs of HQ proteins carry an experimental or TAS/IC code without having been admitted (they count as support here, so goa
can gain a few proteins and hq_explabels 577 labels), and a rule-completed label appears in the annotations only where no
other annotation explains it (so hq_explabels loses those also explained by a computational annotation).
"""
import argparse
import csv

import numpy as np

from retrain.common import Ontology, add_argument, data_root, load, log, read_annotations, write_dataset

EXPERIMENTAL = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'}
CURATOR = {'TAS', 'IC'}
RULE_SOURCE = 'label completion rule (this work)'
PAPER = {'goa': (69_412, 911_848), 'hq_explabels': (87_797, 644_115)}  # proteins, propagated labels


def is_goa(row):
    return row[3] == 'GOA' or row[3].startswith('GOA:')


def support(root):
    """Rows of release/hq_support.tsv by route, or None without the file (with a warning)."""
    path = root / 'release' / 'hq_support.tsv'
    if not path.is_file():
        log(f'WARNING: {path} not found. The set is approximated from the released annotations and will differ slightly '
            "from the paper's (module docstring).")
        return None
    out = {'admitted_goa': [], 'rule_completion': []}
    with open(path, newline='', encoding='utf-8') as handle:
        for r in csv.DictReader(handle, delimiter='\t'):
            source = f"GOA:{r['assigned_by']}" if r['route'] == 'admitted_goa' else RULE_SOURCE
            out[r['route']].append((r['seq_id'], r['go_id'], r['evidence'], source))
    return out


def finish(root, name, ontology, hq, keep, rows, columns=None, exact=False):
    ids = [hq.seq_ids[i] for i in keep]
    labels, dropped = ontology.label_matrix([(h, g) for h, g, _, _ in rows], {h: i for i, h in enumerate(ids)})
    outside = (labels.astype(np.int8) - hq.labels[keep].astype(np.int8)) > 0
    if outside.nnz:
        raise ValueError(f'{name}: {outside.nnz} labels that the same proteins lack in hq')
    found = (len(ids), int(labels.nnz))
    log(f"{name}: {found[0]:,} proteins, {found[1]:,} labels; the paper's: {PAPER[name][0]:,} and {PAPER[name][1]:,} "
        f'({len(dropped)} asserted pairs map to no live term)')
    if exact and found != PAPER[name]:
        raise ValueError(f"{name}: {found} differs from the paper's {PAPER[name]}")
    write_dataset(root, name, ids, [hq.accessions[i].split(',') if hq.accessions[i] else [] for i in keep],
                  [hq.sequences[i] for i in keep], rows, labels, ontology, columns=columns)


def goa(root, ontology, hq, sidecar):
    rows = [r for r in read_annotations(hq.folder) if is_goa(r)]
    admitted = sidecar['admitted_goa'] if sidecar else [r for r in rows if r[2] in EXPERIMENTAL | CURATOR]
    members = {r[0] for r in admitted}
    experimental = {r[0] for r in admitted if r[2] in EXPERIMENTAL}
    keep = [i for i, h in enumerate(hq.seq_ids) if h in members]
    rows = sorted(set(r for r in rows + list(admitted) if r[0] in members))
    log(f'goa: {len(experimental):,} proteins with an experimental code, {len(members) - len(experimental):,} with TAS/IC '
        'only')
    finish(root, 'goa', ontology, hq, keep, rows, dict(experimental=[int(hq.seq_ids[i] in experimental) for i in keep]),
           exact=sidecar is not None)


def hq_explabels(root, ontology, hq, sidecar):
    rows = read_annotations(hq.folder)
    if sidecar:
        rows = [r for r in rows if not is_goa(r) and r[3] != RULE_SOURCE]
        rows += sidecar['admitted_goa'] + sidecar['rule_completion']
    else:
        rows = [r for r in rows if not is_goa(r) or r[2] in EXPERIMENTAL | CURATOR]
    finish(root, 'hq_explabels', ontology, hq, list(range(len(hq.seq_ids))), rows, exact=sidecar is not None)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('dataset', choices=['goa', 'hq_explabels'])
    add_argument(parser)
    a = parser.parse_args()
    root = data_root(a.data)
    {'goa': goa, 'hq_explabels': hq_explabels}[a.dataset](root, Ontology(root), load(root, 'hq'), support(root))


if __name__ == '__main__':
    main()
