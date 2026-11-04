"""Finch and the published methods on their own benchmarks with the test labels' unannotated pairs masked: is a model's
deficit on these experimental-only test labels an artefact of predictions of computationally annotated functions being
counted as errors? The test labels with every current GOA molecular-function annotation of the test proteins masked, as
the paper's "annotations masked" test version masks them (archipelago.evaluate.Benchmark mask: such a pair is neither
positive nor negative); saved predictions only, the paper's evaluator only (Supplementary Table 17).

    python -m retrain.outside.masked fetch      UniProt's secondary accessions (sec_ac.txt, current release) and DPFunc's
                                                entry names mapped to accessions (UniProt ID mapping, REST)
    python -m retrain.outside.masked scan       the test proteins' current accessions (one pass over GOA's GPI) and their
                                                GOA rows (one streamed pass over the GAF)
    python -m retrain.outside.masked prepare METHOD     CPU: the score, truth and mask matrices -> WORK/<method>/arrays_*.npz
    python -m retrain.outside.masked bootstrap METHOD   GPU (archipelago.evaluate on CUDA; --device cpu works, slower):
                                                        -> WORK/<method>/{table,contrasts,masks}.tsv, timing.json
    every stage: --data WORK --outside OUT --replication REP [--goa GOA] [--swissprot FASTA]

WORK (--data) receives everything; OUT/<method> and REP/<method> are each method's outside root (retrain.outside.build,
with Finch combined) and replication root (retrain.replicate); GOA holds goa_uniprot_all.gaf.gz and
goa_uniprot_all.gpi.gz (the paper: the GOA UniProt release of 2026-07-29); FASTA is Swiss-Prot (uniprot_sprot.fasta.gz,
the paper: UniProt 2026_02). figures/si.py reads WORK as <root>/results/outside_masked.

Accessions: DeepGO-SE's records list their entry's accessions, TransFew's one; DPFunc's give only the entry name, mapped
by UniProt's ID mapping (and, for names it no longer knows, the Swiss-Prot entries with the identical sequence in
Swiss-Prot). A secondary accession is replaced by its current primary (sec_ac.txt). A record is matched
when one of its accessions is a current entry of GOA's GPI; the others keep their labels unmasked (counted).
Annotations: UniProtKB rows, aspect F, no NOT qualifier, any evidence; each GO id mapped onto the paper's GO release
(Ontology.map_term) and propagated over the method's own relations (build.RELATIONS), within its vocabulary without the
root (the columns compare.py scores). Evidence: experimental (EXP IDA IPI IMP IGI IEP HTP HDA HMP HGI HEP) and curator
statements (TAS IC) against the rest (IEA IBA ISS ISO ISA ISM IGC RCA NAS ND ...).
Masks: computational = pairs the rest supports (propagated) that experimental or curator support does not, and that are
not positive; all_goa = every pair current GOA supports that is not positive (it also covers experimental annotations
added since the benchmark was frozen). The paper's own positives stay positive.
Metrics: micro-AUPRC and micro-Fmax (archipelago.evaluate.Benchmark), unmasked and with each mask; models, subsets
(all, no_relative) and the paired bootstrap as compare.py (unique test sequences, 10,000 replicates, the same draws for
every model and version, seed means inside each replicate).
"""
import argparse
import csv
import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np

from retrain.common import add_argument, data_root, log
from retrain.outside.build import RELATIONS, ancestors_over
from retrain.outside.compare import CONTRASTS, ROOT_TERM, Bench, read_tsv

# set by main from the command line: the work folder, its inputs, GOA's files, the roots, Swiss-Prot
W = INPUTS = GOA_RAW = OUTSIDE = REPLICATION = SWISSPROT = None
DEVICE = 'cuda'
SEC_AC ='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/sec_ac.txt'
RELDATE = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt'
CURATED = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP', 'TAS', 'IC'}
METHODS = ('deepgose', 'transfew', 'dpfunc')


def bench(method):
    return Bench(method, OUTSIDE / method, REPLICATION / method)


def swissprot_by_sequence(wanted):
    """{sequence SHA-256: {Swiss-Prot accessions}} of the sequences in `wanted` (uppercase SHA-256, as the seq_ids)."""
    import gzip
    import hashlib
    out, accession, parts = {}, None, []

    def flush():
        if accession:
            h = hashlib.sha256(''.join(parts).upper().encode()).hexdigest()
            if h in wanted:
                out.setdefault(h, set()).add(accession)
    with gzip.open(SWISSPROT, 'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                flush()
                accession, parts = line[1:].split('|')[1], []
            else:
                parts.append(line.strip())
    flush()
    return out


def get(url, data=None):
    import urllib.request
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, data=data, timeout=300) as r:
                return r.read()
        except OSError as e:
            log('retry', url, e)
            time.sleep(10 * (attempt + 1))
    raise RuntimeError(f'could not fetch {url}')


def fetch(_):
    import urllib.parse
    INPUTS.mkdir(parents=True, exist_ok=True)
    (INPUTS / 'reldate.txt').write_bytes(get(RELDATE))
    (INPUTS / 'sec_ac.txt').write_bytes(get(SEC_AC))
    names = sorted({r['record'] for r in bench('dpfunc').records})
    job = json.loads(get('https://rest.uniprot.org/idmapping/run', urllib.parse.urlencode(
        {'from': 'UniProtKB_AC-ID', 'to': 'UniProtKB', 'ids': ','.join(names)}).encode()))['jobId']
    for _ in range(120):
        status = json.loads(get(f'https://rest.uniprot.org/idmapping/status/{job}'))
        if status.get('jobStatus') not in ('NEW', 'RUNNING'):
            break
        time.sleep(5)
    table = get(f'https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{job}?format=tsv&fields=accession,id')
    (INPUTS / 'dpfunc_entry_names.tsv').write_bytes(table)
    log('FETCH', (INPUTS / 'reldate.txt').read_text().strip().replace('\n', '; '), 'sec_ac bytes',
        (INPUTS / 'sec_ac.txt').stat().st_size, 'DPFunc names', len(names), 'mapped rows', table.count(b'\n') - 1)


def secondary():
    """{secondary accession: [primary accessions]} from sec_ac.txt."""
    out, started = {}, False
    for line in open(INPUTS / 'sec_ac.txt'):
        if line.startswith('_____'):
            started = True
            continue
        parts = line.split()
        if started and len(parts) == 2:
            out.setdefault(parts[0], []).append(parts[1])
    return out


def stream_grep(source, patterns, out):
    """pigz -dc source | grep -F -w -f patterns > out, with progress."""
    with open(out, 'wb') as handle, open(source, 'rb') as raw:
        unzip = subprocess.Popen(['pigz', '-dc'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        grep = subprocess.Popen(['grep', '-F', '-w', '-f', str(patterns)], stdin=unzip.stdout, stdout=handle)
        unzip.stdout.close()
        total, done, started, shown = os.path.getsize(source), 0, time.monotonic(), time.monotonic()
        for block in iter(lambda: raw.read(64 << 20), b''):
            unzip.stdin.write(block)
            done += len(block)
            if time.monotonic() - shown > 60:
                shown, elapsed = time.monotonic(), time.monotonic() - started
                log('SCAN', source.name, f'{done / total:.1%}', f'eta {elapsed * (total - done) / done / 60:.1f} min')
        unzip.stdin.close()
        if unzip.wait() != 0 or grep.wait() not in (0, 1):
            raise RuntimeError(f'scan of {source} failed')


def scan(_):
    sec = secondary()
    names = {}
    for r in read_tsv(INPUTS / 'dpfunc_entry_names.tsv'):
        names.setdefault(r['From'], []).append(r['Entry'])
    wanted = {}  # (method, record) -> {'listed': [...], 'route': ...}
    for method in METHODS:
        for r in bench(method).records:
            if method == 'dpfunc':
                listed, route = names.get(r['record'], []), 'entry name (UniProt ID mapping)'
            else:
                listed, route = [a for a in r['accession'].split(',') if a], 'accession'
            wanted[method, r['record']] = dict(seq_id=r['seq_id'], listed=listed, route=route)
    # DPFunc names the ID mapping no longer knows: Swiss-Prot entries with the identical sequence
    unknown = {v['seq_id'] for (m, _), v in wanted.items() if m == 'dpfunc' and not v['listed']}
    if unknown and SWISSPROT is None:
        raise SystemExit(f'{len(unknown)} DPFunc entry names are unknown to the ID mapping: give --swissprot')
    by_sequence = swissprot_by_sequence(unknown) if unknown else {}
    candidates = set()
    for v in wanted.values():
        for a in v['listed']:
            candidates.add(a)
            candidates.update(sec.get(a, []))
    for accs in by_sequence.values():
        candidates.update(accs)
    (W / 'scan').mkdir(parents=True, exist_ok=True)
    (W / 'scan' / 'candidates.txt').write_text(''.join(a + '\n' for a in sorted(candidates)))
    stream_grep(GOA_RAW / 'goa_uniprot_all.gpi.gz', W / 'scan' / 'candidates.txt', W / 'scan' / 'gpi_hits.tsv')
    current, swissprot = set(), set()
    for line in open(W / 'scan' / 'gpi_hits.tsv'):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'UniProtKB' and f[1] in candidates:
            current.add(f[1])
            if 'db_subset=Swiss-Prot' in line:
                swissprot.add(f[1])
    matched, counts = {}, {}
    for (method, record), v in wanted.items():
        accs, route = set(), v['route']
        for a in v['listed']:
            if a in current:
                accs.add(a)
            for p in sec.get(a, []):
                if p in current:
                    accs.add(p)
                    route = v['route'] + ', secondary accession'
        if not accs and method == 'dpfunc':
            accs = {a for a in by_sequence.get(v['seq_id'], ()) if a in swissprot}
            route = 'identical sequence (Swiss-Prot)' if accs else route
        matched[method, record] = sorted(accs)
        key = f"{method}: {route if accs else 'not matched'}"
        counts[key] = counts.get(key, 0) + 1
    keep = sorted({a for accs in matched.values() for a in accs})
    (W / 'scan' / 'accessions.txt').write_text(''.join(a + '\n' for a in keep))
    stream_grep(GOA_RAW / 'goa_uniprot_all.gaf.gz', W / 'scan' / 'accessions.txt', W / 'scan' / 'gaf_hits.tsv')
    rows = 0
    keep_set = set(keep)
    with open(W / 'scan' / 'annotations.tsv', 'w') as out:
        out.write('accession\tgo_id\tevidence\tqualifier\n')
        for line in open(W / 'scan' / 'gaf_hits.tsv'):
            f = line.rstrip('\n').split('\t')
            if len(f) >= 15 and f[0] == 'UniProtKB' and f[1] in keep_set and f[8] == 'F' and not f[3].startswith('NOT'):
                out.write(f'{f[1]}\t{f[4]}\t{f[6]}\t{f[3]}\n')
                rows += 1
    (W / 'scan' / 'matched.json').write_text(json.dumps({f'{m}\t{r}': a for (m, r), a in matched.items()}) + '\n')
    record = dict(candidates=len(candidates), current_entries=len(current), accessions_scanned=len(keep),
                  mf_annotation_rows=rows, records=counts, gaf=str(GOA_RAW / 'goa_uniprot_all.gaf.gz'),
                  uniprot_release=(INPUTS / 'reldate.txt').read_text().strip())
    (W / 'scan' / 'scan.json').write_text(json.dumps(record, indent=1) + '\n')
    log('SCAN', json.dumps(record))


def prepare(a):
    """CPU: the paper's propagated prediction and truth matrices (compare.py's), the two masks, and the rows' sequences and
    relatives -> <W>/<method>/arrays_<test set>.npz."""
    from retrain.common import Ontology
    from retrain.outside import compare
    method = a.method
    b = bench(method)
    views = getattr(compare, method)(b, with_counts=False)
    matched = {k.split('\t', 1)[1]: v for k, v in json.loads((W / 'scan' / 'matched.json').read_text()).items()
               if k.split('\t', 1)[0] == method}
    support = {}  # accession -> {'curated': set, 'all': set}
    for r in read_tsv(W / 'scan' / 'annotations.tsv'):
        s = support.setdefault(r['accession'], {'curated': set(), 'all': set()})
        s['all'].add(r['go_id'])
        if r['evidence'] in CURATED:
            s['curated'].add(r['go_id'])
    ontology = Ontology(b.path)
    ancestors = ancestors_over(b.path / 'datasets' / 'ontology' / 'go-basic.obo', RELATIONS[method])
    out = W / method
    out.mkdir(parents=True, exist_ok=True)
    for test_set, (rows, _, (truth, scores, terms), _) in views.items():
        column = {t: j for j, t in enumerate(terms)}

        def closed(go_ids):
            cols = set()
            for g in go_ids:
                t = ontology.map_term(g)
                if t is not None:
                    cols.update(column[x] for x in ancestors(t) if x in column)
            return sorted(cols)
        curated, everything = np.zeros(truth.shape, bool), np.zeros(truth.shape, bool)
        for i, r in enumerate(rows):
            accs = matched.get(r['record'], [])
            curated[i, closed(set().union(set(), *(support.get(x, {}).get('curated', set()) for x in accs)))] = True
            everything[i, closed(set().union(set(), *(support.get(x, {}).get('all', set()) for x in accs)))] = True
        keys = sorted(scores)
        np.savez(out / f'arrays_{test_set}.npz', truth=truth, computational=everything & ~curated & ~truth,
                 all_goa=everything & ~truth, seq=np.array([r['seq_id'] for r in rows]),
                 record=np.array([r['record'] for r in rows]),
                 relative=np.array([r['relative'] for r in rows]),
                 matched=np.array([bool(matched.get(r['record'])) for r in rows]),
                 keys=np.array([f'{m}|{s}' for m, s in keys]), **{f'scores_{i}': scores[k] for i, k in enumerate(keys)})
        log(method, test_set, 'rows', len(rows), 'not matched to GOA', int(sum(not matched.get(r['record']) for r in rows)),
            'masked pairs: computational', int((everything & ~curated & ~truth).sum()), 'all GOA', int((everything & ~truth).sum()))


def bootstrap(a):
    """GPU: micro-AUPRC and micro-Fmax (archipelago.evaluate.Benchmark on CUDA), unmasked and with each mask, 10,000
    paired replicates; -> <W>/<method>/{table,contrasts,masks}.tsv and timing.json."""
    from archipelago.evaluate import Benchmark, cluster_draws
    started = time.monotonic()
    method, out = a.method, W / a.method
    table, contrasts, masks_out, timing = [], [], [], {}
    for path in sorted(out.glob('arrays_*.npz')):
        test_set = path.stem[len('arrays_'):]
        z = np.load(path)
        truth, seq, relative = z['truth'], z['seq'], z['relative']
        keys = [tuple(k.split('|')) for k in z['keys']]
        scores = {k: z[f'scores_{i}'] for i, k in enumerate(keys)}
        models = list(dict.fromkeys(m for m, _ in keys))
        versions = {'unmasked': None, 'computational': z['computational'], 'all_goa': z['all_goa']}
        for subset, mask in (('all', np.ones(len(seq), bool)), ('no_relative', relative == 'none')):
            rows_in = np.flatnonzero(mask)
            index, draws = cluster_draws(seq[rows_in])
            labelled = rows_in[truth[rows_in].any(1)]
            present, relabel = np.unique(index[np.isin(rows_in, labelled)], return_inverse=True)
            for version, m in versions.items():
                t0 = time.monotonic()
                bench = Benchmark(truth[labelled], clusters=relabel, draws=draws[:, present], device=DEVICE,
                                  mask=None if m is None else m[labelled])
                masked_pairs = 0 if m is None else int(m[labelled].sum())
                masks_out.append(dict(test_set=test_set, subset=subset, version=version, proteins=len(labelled),
                                      proteins_not_matched_to_GOA=int((~z['matched'][labelled]).sum()),
                                      masked_pairs=masked_pairs, positive_pairs=int(truth[labelled].sum())))
                stats = {}
                for model in models:
                    ks = [k for k in keys if k[0] == model]
                    res = [bench.score(scores[k][labelled]) for k in ks]
                    au = np.array([r['metrics']['auprc'] for r in res])
                    fm = np.array([r['metrics']['fmax'] for r in res])
                    stats[model] = dict(micro_AUPRC=(au, np.mean([r['replicates'][:, 0] for r in res], axis=0)),
                                        micro_Fmax=(fm, np.mean([r['replicates'][:, 1] for r in res], axis=0)))
                    sd = lambda x: float(x.std(ddof=1)) if len(x) > 1 else 0.0  # noqa: E731
                    table.append(dict(test_set=test_set, subset=subset, version=version, model=model, runs=len(ks),
                                      proteins=len(labelled), sequences=len(present), masked_pairs=masked_pairs,
                                      micro_AUPRC_mean=au.mean(), micro_AUPRC_sd=sd(au), micro_Fmax_mean=fm.mean(),
                                      micro_Fmax_sd=sd(fm)))
                    log(method, table[-1])
                for x, y in CONTRASTS:
                    if x not in stats or y not in stats:
                        continue
                    for metric in ('micro_AUPRC', 'micro_Fmax'):
                        (px, rx), (py, ry) = stats[x][metric], stats[y][metric]
                        low, high = np.percentile(rx - ry, [2.5, 97.5])
                        contrasts.append(dict(test_set=test_set, subset=subset, version=version, contrast=f'{x} - {y}',
                                              metric=metric, difference=px.mean() - py.mean(), ci_low=low, ci_high=high,
                                              n_units=len(present), replicates=len(draws)))
                timing[f'{test_set}/{subset}/{version}'] = round(time.monotonic() - t0, 1)
    timing['total_seconds'] = round(time.monotonic() - started, 1)
    for name, data in (('table', table), ('contrasts', contrasts), ('masks', masks_out)):
        with open(out / f'{name}.tsv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(data[0]), delimiter='\t')
            w.writeheader()
            w.writerows(data)
    (out / 'timing.json').write_text(json.dumps(timing, indent=1) + '\n')
    log(method, 'wrote', out, json.dumps(timing))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('stage', choices=['fetch', 'scan', 'prepare', 'bootstrap'])
    parser.add_argument('method', nargs='?', choices=METHODS)
    parser.add_argument('--outside', type=Path, required=True, help="the outside roots' parent (OUT/<method>)")
    parser.add_argument('--replication', type=Path, required=True, help="the replication roots' parent (REP/<method>)")
    parser.add_argument('--goa', type=Path, help='the folder of goa_uniprot_all.gaf.gz and goa_uniprot_all.gpi.gz (scan)')
    parser.add_argument('--swissprot', type=Path, help='uniprot_sprot.fasta.gz (scan, for DPFunc names the ID mapping lacks)')
    parser.add_argument('--device', default='cuda', help='bootstrap: cuda (the paper) or cpu')
    add_argument(parser)
    a = parser.parse_args()
    global W, INPUTS, GOA_RAW, OUTSIDE, REPLICATION, SWISSPROT, DEVICE
    W = data_root(a.data)
    INPUTS = W / 'inputs'
    GOA_RAW, OUTSIDE, REPLICATION, SWISSPROT, DEVICE = a.goa, a.outside, a.replication, a.swissprot, a.device
    if a.stage == 'scan' and GOA_RAW is None:
        raise SystemExit('scan needs --goa')
    if a.stage in ('prepare', 'bootstrap') and a.method is None:
        raise SystemExit(f'{a.stage} needs a method')
    dict(fetch=fetch, scan=scan, prepare=prepare, bootstrap=bootstrap)[a.stage](a)


if __name__ == '__main__':
    main()
