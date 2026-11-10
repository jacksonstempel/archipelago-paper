"""The sample of the second review: 200 accepted and 200 not-accepted test assertions, the first review's final
decisions. The 200 of the first draw (first_sample.tsv) plus 100 more of each from the rest (seed 20261031); all 400 shuffled
and split into 16 packets of 25, built like the original review files (the original PROMPT.md unchanged; each protein's
sequence, identifiers and assertion records as the original packets gave them, JSON indented). The key is written
separately.

    python make_packets.py REVIEW_DIR RESOLVED_JSON first_sample.tsv PACKET_DIR KEY_TSV

The second review was run on packets 1 to 4 (100 assertions; packets/ holds them as sent, results/ what came back). The
inputs REVIEW_DIR and RESOLVED_JSON are the review's own files (make_first_sample.py), not part of the release.
"""
import csv
import json
import random
import sys
from pathlib import Path

review_dir, resolved, first_key, packet_dir, key = (Path(a) for a in sys.argv[1:6])
SEED, MORE_PER_CLASS, PACKETS = 20261031, 100, 16

decisions = {a['assertion_id']: a['decision'] for a in json.load(open(resolved, encoding='utf-8'))}
proteins, assertion = {}, {}
for f in sorted(review_dir.glob('PRO_TEST_REVIEW_*_of_41.md')):
    batch = json.loads(f.read_text(encoding='utf-8').split('```json\n', 1)[1].rsplit('\n```', 1)[0])
    for p in batch['proteins']:
        proteins[p['sequence_sha256']] = {k: v for k, v in p.items() if k != 'assertions'}
        for a in p['assertions']:
            assertion[a['assertion_id']] = (p['sequence_sha256'], a)
if set(assertion) != set(decisions):
    raise ValueError('the packets and the final decisions list different assertions')

first = [(r['assertion_id'], r['original_decision']) for r in csv.DictReader(open(first_key, encoding='utf-8'), delimiter='\t')]
taken = {i for i, _ in first}
rng = random.Random(SEED)
sample = list(first)
for decision in ('accept', 'do_not_accept'):
    pool = sorted(i for i, d in decisions.items() if d == decision and i not in taken)
    sample += [(i, decision) for i in rng.sample(pool, MORE_PER_CLASS)]
rng.shuffle(sample)

prompt = (review_dir / 'PROMPT.md').read_text(encoding='utf-8')
packet_dir.mkdir(parents=True)
rows = []
for k in range(PACKETS):
    part = sample[k * len(sample) // PACKETS:(k + 1) * len(sample) // PACKETS]
    grouped = {}
    for i, decision in part:
        seq, a = assertion[i]
        grouped.setdefault(seq, []).append(a)
        rows.append(dict(packet=k + 1, assertion_id=i, original_decision=decision, in_first_check=int(i in taken)))
    data = dict(batch=k + 1, source_release='GOA 2026-07-29',
                proteins=[dict(proteins[seq], assertions=a) for seq, a in grouped.items()])
    name = f'AGREEMENT_RESULTS_{k + 1:02d}.json'
    text = (f'# Protein function review — batch {k + 1} of {PACKETS}\n\n{len(grouped)} proteins; {len(part)} '
            f'assertions. Name your output `{name}`.\n\n{prompt.rstrip()}\n\n\n## Protein and assertion data\n\n```json\n'
            + json.dumps(data, ensure_ascii=False, indent=1) + '\n```\n')
    (packet_dir / f'AGREEMENT_REVIEW_{k + 1:02d}_of_{PACKETS}.md').write_text(text, encoding='utf-8')
    print(f'packet {k + 1}: {len(grouped)} proteins, {len(part)} assertions, '
          f'{sum(d == "accept" for _, d in part)} originally accepted')
with open(key, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader()
    w.writerows(rows)
