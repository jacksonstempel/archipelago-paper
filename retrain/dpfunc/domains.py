"""DPFunc's domain input: InterPro entries found by InterProScan 5.78-109.0, counted per protein.

    python -m retrain.dpfunc.domains scan --interproscan DIR [--cpus 4] [--part i --parts n]
    python -m retrain.dpfunc.domains merge

scan   InterProScan (all 18 member analyses of release 5.78-109.0, with the precalculated-match lookup) on batches of
       1,024 proteins of the datasets (sorted by length); batch k belongs to part k mod n. Writes
       inputs/dpfunc/domains/batches/<k>/counts.json: per protein, how many match rows name each integrated InterPro entry
       (repeated entries counted, as the released DPFunc preprocessing does). Finished batches are skipped. About 90 s
       per 256 proteins on 4 CPUs in the paper's runs.
merge  one table for every protein: inputs/dpfunc/domains/counts.json.gz.

Every match row is checked against its sequence (MD5 and length), and the run must state release 5.78-109.0 with the
analyses below; another InterProScan release gives other entries. InterProScan is Apache 2.0 licensed
(github.com/ebi-pf-team/interproscan); it needs Java 11 and its data files.
"""
import argparse
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import time
from pathlib import Path

from retrain.common import add_argument, data_root, load, log
from retrain.dpfunc.common import inputs

DATASETS = ('original_dpfunc', 'archipelago', 'development', 'test')
BATCH = 1024
ANALYSES = ('[AntiFam-8.0,CDD-3.21,Coils-2.2.1,FunFam-4.3.0,Gene3D-4.3.0,Hamap-2026_01,MobiDBLite-4.0,NCBIfam-19.0,'
            'PANTHER-19.0,Pfam-38.2,PIRSF-3.10,PIRSR-2025_05,PRINTS-42.0,ProSitePatterns-2026_01,ProSiteProfiles-2026_01,'
            'SFLD-4,SMART-9.0,SUPERFAMILY-1.75]')


def proteins(root, names=DATASETS):
    sequences = {}
    for name in names:
        ds = load(root, name, labels=False)
        sequences.update(zip(ds.seq_ids, ds.sequences))
    return sequences


def batches(sequences):
    ordered = sorted(sequences, key=lambda h: (len(sequences[h]), h))
    return [ordered[a:a + BATCH] for a in range(0, len(ordered), BATCH)]


def analysis_versions(text):
    found = [v for v in re.findall(r'Running the following analyses:\s*\n(\[[^\n]+\])', text) if not v.startswith('[job')]
    if len(set(found)) != 1 or 'InterProScan-5.78-109.0' not in text:
        raise RuntimeError('the scan log does not state one set of InterProScan 5.78-109.0 analyses')
    return found[0]


def scan_counts(sequences, matches):
    """{seq_id: {InterPro entry: number of match rows}}; every row is checked by MD5 and length."""
    counts = {h: {} for h in sequences}
    for line in Path(matches).read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        sequence = sequences[fields[0]]
        if fields[1].lower() != hashlib.md5(sequence.encode()).hexdigest() or int(fields[2]) != len(sequence):
            raise RuntimeError(f'{matches}: a match row does not fit its sequence')
        if len(fields) > 11 and fields[11].startswith('IPR'):
            counts[fields[0]][fields[11]] = counts[fields[0]].get(fields[11], 0) + 1
    return counts


def scan(root, arguments):
    sequences = proteins(root, arguments.datasets)
    folder = inputs(root) / 'domains' / 'batches'
    mine = [(k, b) for k, b in enumerate(batches(sequences)) if k % arguments.parts == arguments.part]
    todo = [(k, b) for k, b in mine if not (folder / f'{k:05d}' / 'counts.json').exists()]
    log(f'InterProScan: {len(todo)} of {len(mine)} batches to scan')
    interproscan = Path(arguments.interproscan).resolve()  # absolute: the scan runs inside this folder
    started = time.monotonic()
    for i, (k, batch) in enumerate(todo, 1):
        work = folder / f'{k:05d}'
        shutil.rmtree(work, ignore_errors=True)
        (work / 'tmp').mkdir(parents=True)
        (work / 'input.fasta').write_text(''.join(f'>{h}\n{sequences[h]}\n' for h in batch))
        with open(work / 'scan.log', 'w') as out:
            subprocess.run([str(interproscan / 'interproscan.sh'), '-i', work / 'input.fasta', '-f', 'tsv',
                            '-iprlookup', '-o', work / 'matches.tsv', '-cpu', str(arguments.cpus), '-T', work / 'tmp'],
                           cwd=interproscan, stdout=out, stderr=subprocess.STDOUT, check=True)
        text = (work / 'scan.log').read_text()
        if '100% done:  InterProScan analyses completed' not in text or analysis_versions(text) != ANALYSES:
            raise RuntimeError(f'{work}: incomplete scan, or analyses other than release 5.78-109.0')
        counts = scan_counts({h: sequences[h] for h in batch}, work / 'matches.tsv')
        shutil.rmtree(work / 'tmp')
        (work / 'counts.json').write_text(json.dumps(counts))
        rate = i / (time.monotonic() - started)
        log(f'batch {k}: {len(batch)} proteins; {i}/{len(todo)} batches, ETA {(len(todo) - i) / rate / 3600:.1f} h')


def merge(root, arguments):
    sequences = proteins(root, arguments.datasets)
    folder = inputs(root) / 'domains'
    counts = {}
    for k, batch in enumerate(batches(sequences)):
        path = folder / 'batches' / f'{k:05d}' / 'counts.json'
        if not path.exists():
            raise SystemExit(f'batch {k} has not been scanned; run the scan stage')
        found = json.loads(path.read_text())
        if set(found) != set(batch):
            raise RuntimeError(f'{path} does not cover its batch (were the datasets changed after the scan?)')
        counts.update(found)
    with gzip.open(folder / 'counts.json.gz', 'wt') as handle:
        json.dump(dict(sorted(counts.items())), handle)
    log(f'domains: {len(counts):,} proteins, {sum(not c for c in counts.values()):,} without an InterPro entry')


def read_counts(root):
    with gzip.open(inputs(root) / 'domains' / 'counts.json.gz', 'rt') as handle:
        return json.load(handle)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('stage', choices=['scan', 'merge'])
    parser.add_argument('--datasets', nargs='+', default=list(DATASETS))
    parser.add_argument('--interproscan', help='scan: the InterProScan 5.78-109.0 folder (holding interproscan.sh)')
    parser.add_argument('--cpus', type=int, default=4)
    parser.add_argument('--part', type=int, default=0)
    parser.add_argument('--parts', type=int, default=1)
    add_argument(parser)
    arguments = parser.parse_args()
    if arguments.stage == 'scan' and not arguments.interproscan:
        parser.error('scan needs --interproscan')
    globals()[arguments.stage](data_root(arguments.data), arguments)


if __name__ == '__main__':
    main()
