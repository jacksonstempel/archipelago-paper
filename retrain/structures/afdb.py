"""Structure bank, part 1: an AlphaFold DB model of each protein's exact sequence.

    python -m retrain.structures.afdb requests [--datasets original_dpfunc archipelago development test]
    python -m retrain.structures.afdb fetch --chosen [FILE]      (the paper's bank; default inputs/structures/chosen.tsv)
    python -m retrain.structures.afdb search [--limit N]; python -m retrain.structures.afdb fetch [--limit N]   (free search)

requests  every distinct protein of the datasets -> requests.tsv, requests.fasta (in <data root>/inputs/structures/).
          The defaults cover every set DPFunc is trained or scored on (goa and hq are subsets of archipelago).
fetch --chosen [FILE]
          the paper's choice: FILE (released with the data; columns seq_id, source = afdb or colabfold, entry, version)
          names for every protein either the AlphaFold DB model the paper used (e.g. AF-P12345-F1, 6) or ColabFold. Exactly
          those models are downloaded and verified; a protein missing from FILE, or a listed model that cannot be
          downloaded or does not match the protein, stops the stage.
search    for a free search instead: AlphaFold DB's sequence search (by the MD5 of the modeled residues: the whole
          sequence, or residues 1-1,400 above 2,700 residues) -> afdb/search.jsonl. A sequence the search answers with
          HTTP 500 (it fails reproducibly for some sequences) is looked up through every UniProt accession of the identical
          sequence (UniParc, by CRC64) and AlphaFold DB's models of those accessions. Resumable.
fetch     (free search) downloads one verified model per protein that has one: model versions up to 6 (the release the
          paper used), the newest first, then one of the protein's own UniProt accessions, then the entry id.

Both fetch modes verify each model (residues identical to the modeled residues, complete N/CA/C) and write afdb/files/,
afdb/selected.tsv and, for every other protein, needs_colabfold.tsv (folded by retrain.structures.colabfold). Resumable.

How a free search differs from the paper's bank, which also took models from local archives: it does not find the
version-4 F1 fragments of long human proteins that the paper's bank kept where no version-6 model exists, and the
paper's search also sent a control sample of HTTP-404 answers through the UniParc route (to check that 404 means no
model). Use --chosen to reproduce the paper.

Only Google DeepMind's AlphaFold models are used (AF-<UniProt accession>-F<n>); AlphaFold DB also hosts models from other
providers, which are ignored. AlphaFold DB data are CC BY 4.0 (cite Fleming et al. 2025 and Jumper et al. 2021).
"""
import argparse
import hashlib
import json
import re
import threading
import time
import urllib.parse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from retrain.common import add_argument, data_root, load, log, read_fasta
from retrain.structures.common import (Progress, bank, modeled_raw, read_jsonl, read_structure, read_tsv, target_sequence,
                                       write_bytes, write_tsv)

API = 'https://alphafold.ebi.ac.uk/api/sequence/summary?id={md5}&type=md5'
PREDICTION = 'https://alphafold.ebi.ac.uk/api/prediction/{accession}'
FILES = 'https://alphafold.ebi.ac.uk/files/{entry}-model_v{version}.pdb'
UNIPARC = 'https://rest.uniprot.org/uniparc/stream?format=tsv&fields=upi,accession,checksum,length&query='
DEEPMIND = re.compile(r'^AF-[A-Z0-9]+(-\d+)?-F\d+$')
DATASETS = ('original_dpfunc', 'archipelago', 'development', 'test')
MAX_VERSION = 6  # the AlphaFold DB release the paper used; a later release would change the inputs
WORKERS = 8  # concurrent API requests
FILE_RATE = 2.0  # files per second to start with (the file host answers HTTP 429 above about 2 per second)


def load_requests(root):
    rows = read_tsv(bank(root) / 'requests.tsv')
    ids, sequences = read_fasta(bank(root) / 'requests.fasta')
    sequence = dict(zip(ids, sequences))
    for r in rows:
        r['sequence'] = sequence[r['seq_id']]
        r['own'] = set(filter(None, r['accessions'].split(',')))
    return rows


def requests(root, arguments):
    proteins = {}
    for name in arguments.datasets:
        ds = load(root, name, labels=False)
        for h, accessions, sequence in zip(ds.seq_ids, ds.accessions, ds.sequences):
            p = proteins.setdefault(h, dict(sequence=sequence.upper(), sets=set(), accessions=set()))
            p['sets'].add(name)
            p['accessions'].update(filter(None, accessions.split(',')))
    folder = bank(root)
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / 'requests.fasta', 'w') as handle:
        handle.writelines(f'>{h}\n{p["sequence"]}\n' for h, p in sorted(proteins.items()))
    write_tsv(folder / 'requests.tsv', [dict(seq_id=h, length=len(p['sequence']), sets=','.join(sorted(p['sets'])),
                                             accessions=','.join(sorted(p['accessions']))) for h, p in sorted(proteins.items())],
              ['seq_id', 'length', 'sets', 'accessions'])
    log(f'requests: {len(proteins):,} distinct proteins from {", ".join(arguments.datasets)}')


def api_get(session, url, attempts=8, server_errors=3):
    """(status, JSON or None). 404 means no model. 429 and network errors are retried with backoff; a 500 repeated
    `server_errors` times is returned as 500."""
    import requests as http
    errors = 0
    for attempt in range(attempts):
        try:
            response = session.get(url, timeout=(20, 60))
            if response.status_code == 200:
                return 200, response.json()
            if response.status_code == 404:
                return 404, None
            if response.status_code == 500:
                errors += 1
                if errors >= server_errors:
                    return 500, None
            delay = float(response.headers.get('Retry-After', 2 ** attempt)) if response.status_code == 429 else 2 ** attempt
        except (http.exceptions.RequestException, ValueError):
            delay = 2 ** attempt
        time.sleep(min(delay, 120))
    raise RuntimeError(f'AlphaFold DB API failed {attempts} times: {url}')


def entries_of(answer):
    out = []
    for s in (answer or {}).get('structures', []):
        s = s['summary']
        if s.get('provider') != 'AlphaFold DB' or not DEEPMIND.match(s['model_identifier']):
            continue
        accession = [e['identifier'] for e in s.get('entities', []) if e.get('identifier_category') == 'UNIPROT']
        out.append(dict(entry=s['model_identifier'], version=int(s['model_url'].rsplit('_v', 1)[1].split('.')[0]),
                        accession=accession[0] if accession else ''))
    return out


def crc64(sequence):
    """UniProt's CRC64 checksum (ISO polynomial), as UniParc stores it."""
    table = []
    for i in range(256):
        c = i
        for _ in range(8):
            c = (c >> 1) ^ 0xD800000000000000 if c & 1 else c >> 1
        table.append(c)
    crc = 0
    for ch in sequence.encode():
        crc = table[(crc ^ ch) & 0xff] ^ (crc >> 8)
    return '%016X' % crc


def by_accessions(session, sequence):
    """AlphaFold DeepMind models of the UniProt accessions whose sequence is exactly `sequence` (UniParc by CRC64), kept
    when their modeled residues equal the protein's modeled residues."""
    import requests as http
    raw, _ = modeled_raw(sequence)
    text = None
    query = urllib.parse.quote(f'checksum:{crc64(sequence)}')
    for attempt in range(8):
        try:
            response = session.get(UNIPARC + query, timeout=(20, 300))
            if response.status_code == 200:
                text = response.text
                break
        except http.exceptions.RequestException:
            pass
        time.sleep(min(2 ** attempt, 60))
    if text is None:
        raise RuntimeError('UniParc query failed')
    accessions = set()
    for line in text.splitlines()[1:]:
        _, found, _, length = line.split('\t')
        if int(length) == len(sequence):
            accessions |= {a.split('.')[0] for a in found.split('; ') if a}
    out = []
    for accession in sorted(accessions):
        _, models = api_get(session, PREDICTION.format(accession=accession))
        for m in models or []:
            whole = m.get('uniprotSequence') or m.get('sequence') or ''
            start = int(m.get('uniprotStart') or m.get('sequenceStart') or 1)
            end = int(m.get('uniprotEnd') or m.get('sequenceEnd') or len(whole))
            modeled = m.get('sequence') or ''
            if len(modeled) != end - start + 1:
                modeled = whole[start - 1:end]
            if modeled == raw and (len(raw) == len(sequence) or (whole == sequence and start == 1)) \
                    and DEEPMIND.match(m['entryId']):
                out.append(dict(entry=m['entryId'], version=int(m['latestVersion']),
                                accession=m.get('uniprotAccession', accession)))
    return out


def search(root, arguments):
    import requests as http
    out = bank(root) / 'afdb' / 'search.jsonl'
    out.parent.mkdir(parents=True, exist_ok=True)
    done = {x['seq_id'] for x in read_jsonl(out) if x['status'] != 'error'}
    todo = [r for r in load_requests(root) if r['seq_id'] not in done][:arguments.limit]
    log(f'search: {len(todo):,} proteins to search, {len(done):,} already searched')
    session = http.Session()
    session.mount('https://', http.adapters.HTTPAdapter(pool_connections=WORKERS, pool_maxsize=WORKERS))

    def one(r):
        raw, coverage = modeled_raw(r['sequence'])
        md5 = hashlib.md5(raw.encode()).hexdigest()
        try:
            status, answer = api_get(session, API.format(md5=md5))
            if answer is not None and answer.get('entry', {}).get('sequence') not in (None, raw):
                raise RuntimeError('the API answered for another sequence')
            entries = entries_of(answer) if status != 500 else by_accessions(session, r['sequence'])
        except RuntimeError as error:  # recorded; the next run of this stage tries again
            return dict(seq_id=r['seq_id'], query=coverage, status='error', error=str(error), entries=[])
        return dict(seq_id=r['seq_id'], query=coverage, status=status, entries=entries)

    lock, counts, progress = threading.Lock(), Counter(), Progress('search', len(todo))
    with open(out, 'a') as handle, ThreadPoolExecutor(WORKERS) as pool:
        for future in as_completed([pool.submit(one, r) for r in todo]):
            result = future.result()
            with lock:
                handle.write(json.dumps(result) + '\n')
                handle.flush()
                counts['with_model' if result['entries'] else f"none_{result['status']}"] += 1
                progress.step(**counts)
    if counts.get('none_error'):
        raise SystemExit(f"{counts['none_error']} searches failed; run the search stage again")


class Limiter:
    """Pacing of file downloads: start at FILE_RATE per second, +0.25 after 20 successes up to 4, halve on HTTP 429."""

    def __init__(self):
        self.rate, self.ceiling, self.next, self.pause, self.streak = FILE_RATE, 4.0, 0.0, 0.0, 0
        self.condition = threading.Condition()

    def wait(self):
        with self.condition:
            while (delay := max(self.next, self.pause) - time.monotonic()) > 0:
                self.condition.wait(delay)
            self.next = time.monotonic() + 1 / self.rate

    def success(self):
        with self.condition:
            self.streak += 1
            if self.streak >= 20 and self.rate < self.ceiling:
                self.streak, self.rate = 0, min(self.ceiling, self.rate + 0.25)

    def throttled(self, retry_after):
        with self.condition:
            self.streak = 0
            self.ceiling = max(0.5, self.rate * 0.9)
            self.rate = max(0.5, self.rate / 2)
            self.pause = max(self.pause, time.monotonic() + retry_after)


def download(session, limiter, url):
    import requests as http
    for attempt in range(10):
        limiter.wait()
        try:
            response = session.get(url, timeout=(20, 120))
        except http.exceptions.RequestException:
            time.sleep(min(2 ** attempt, 60))
            continue
        if response.status_code == 200:
            limiter.success()
            return response.content
        if response.status_code == 404:
            return None
        if response.status_code == 429:
            limiter.throttled(float(response.headers.get('Retry-After', 60)))
        else:
            time.sleep(min(2 ** attempt, 60))
    raise RuntimeError(f'download failed 10 times: {url}')


def fetch(root, arguments):
    import gzip

    import requests as http
    folder = bank(root) / 'afdb'
    folder.mkdir(parents=True, exist_ok=True)
    rows = {r['seq_id']: r for r in load_requests(root)}
    chosen = None if arguments.chosen is None else Path(arguments.chosen or bank(root) / 'chosen.tsv')
    candidates = chosen_candidates(chosen, rows) if chosen else searched_candidates(folder, rows)
    done = {r['seq_id'] for r in read_tsv(folder / 'fetched.tsv')} if (folder / 'fetched.tsv').exists() else set()
    todo = [dict(seq_id=h, entries=e) for h, e in candidates.items() if e and h not in done][:arguments.limit]
    log(f'fetch: {len(todo):,} proteins with a model to download, {len(done):,} done')
    columns = ['seq_id', 'status', 'entry', 'version', 'coverage', 'path', 'rejected']
    limiter, local, lock, progress = Limiter(), threading.local(), threading.Lock(), Progress('fetch', len(todo))

    def one(x):
        r = rows[x['seq_id']]
        raw, coverage = modeled_raw(r['sequence'])
        rejected = []
        for e in sorted(x['entries'], key=lambda e: (-e['version'], e['accession'] not in r['own'], e['entry'])):
            path = folder / 'files' / e['entry'][3:5] / f"{e['entry']}-model_v{e['version']}.pdb.gz"
            if path.exists():
                payload = path.read_bytes()
            else:
                if not hasattr(local, 'session'):
                    local.session = http.Session()
                payload = download(local.session, limiter, FILES.format(**e))
                if payload is None:
                    rejected.append(f"{e['entry']}:missing")
                    continue
                payload = gzip.compress(payload, 6)
            try:
                sequence, _, _ = read_structure(payload)
            except ValueError as error:
                rejected.append(f"{e['entry']}:{error}")
                continue
            if sequence != raw:
                rejected.append(f"{e['entry']}:sequence_differs")
                continue
            if not path.exists():
                write_bytes(path, payload)
            return dict(seq_id=x['seq_id'], status='found', entry=e['entry'], version=e['version'], coverage=coverage,
                        path=str(path.relative_to(bank(root))), rejected=';'.join(rejected))
        return dict(seq_id=x['seq_id'], status='no_verified_model', rejected=';'.join(rejected))

    new = not (folder / 'fetched.tsv').exists()
    with open(folder / 'fetched.tsv', 'a') as handle, ThreadPoolExecutor(8) as pool:
        if new:
            handle.write('\t'.join(columns) + '\n')
        for future in as_completed([pool.submit(one, x) for x in todo]):
            result = future.result()
            with lock:
                handle.write('\t'.join(str(result.get(c, '')) for c in columns) + '\n')
                handle.flush()
                progress.step(rate=limiter.rate)
    if arguments.limit:
        return
    fetched = {r['seq_id']: r for r in read_tsv(folder / 'fetched.tsv')}
    if chosen:
        failed = [h for h, e in candidates.items() if e and (fetched[h]['status'] != 'found' or (
            fetched[h]['entry'], int(fetched[h]['version'])) != (e[0]['entry'], e[0]['version']))]
        if failed:
            raise SystemExit(f'{len(failed):,} models of the chosen list could not be verified (first {failed[0]}: '
                             f"{fetched[failed[0]]['rejected']}); see afdb/fetched.tsv")
    selected = [fetched[h] for h in rows if fetched.get(h, {}).get('status') == 'found']
    needs = [dict(seq_id=h, length=len(r['sequence']), target=hashlib.sha256(target_sequence(r['sequence']).encode()).hexdigest(),
                  sets=r['sets']) for h, r in rows.items() if fetched.get(h, {}).get('status') != 'found']
    write_tsv(folder / 'selected.tsv', selected, columns[:-1])
    write_tsv(bank(root) / 'needs_colabfold.tsv', needs, ['seq_id', 'length', 'target', 'sets'])
    log(f'fetch: {len(selected):,} proteins have an AlphaFold DB model, {len(needs):,} need ColabFold')


def chosen_candidates(path, rows):
    """{seq_id: [the one listed model] or [] for ColabFold} from the paper's list; every requested protein must be in it."""
    listed = {r['seq_id']: r for r in read_tsv(path)}
    missing = [h for h in rows if h not in listed]
    if missing:
        raise SystemExit(f'{len(missing):,} requested proteins are not in {path} (first {missing[0]})')
    out = {}
    for h in rows:
        r = listed[h]
        if r['source'] not in ('afdb', 'colabfold'):
            raise ValueError(f"{path}: source {r['source']!r} for {h}")
        out[h] = [dict(entry=r['entry'], version=int(r['version']), accession='')] if r['source'] == 'afdb' else []
    return out


def searched_candidates(folder, rows):
    """{seq_id: DeepMind models of versions up to MAX_VERSION} from the free search."""
    searched = {x['seq_id']: x for x in read_jsonl(folder / 'search.jsonl') if x['status'] != 'error'}
    unsearched = [h for h in rows if h not in searched]
    if unsearched:
        raise SystemExit(f'{len(unsearched):,} proteins were not searched yet; run the search stage first')
    return {h: [e for e in searched[h]['entries'] if e['version'] <= MAX_VERSION] for h in rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('stage', choices=['requests', 'search', 'fetch'])
    parser.add_argument('--datasets', nargs='+', default=list(DATASETS), help='requests: the datasets to cover')
    parser.add_argument('--chosen', nargs='?', const='', metavar='FILE',
                        help="fetch: the paper's list of models (seq_id, source, entry, version); "
                             "default inputs/structures/chosen.tsv")
    parser.add_argument('--limit', type=int, help='search, fetch: only this many proteins (a quick first check)')
    add_argument(parser)
    arguments = parser.parse_args()
    globals()[arguments.stage](data_root(arguments.data), arguments)


if __name__ == '__main__':
    main()
