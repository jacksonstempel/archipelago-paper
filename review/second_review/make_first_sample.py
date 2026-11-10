"""Agreement check of the test-label review: 100 accepted and 100 not-accepted test assertions (final decisions,
resolved_reviews.json), sampled with a fixed seed, shuffled and split into four packets of 50. Each packet is built like
the original review files: the same prompt (PROMPT.md, unchanged) and, for each protein, the sequence, identifiers and
the assertion records exactly as the original packets gave them (JSON indented, one field per line). The key (original
decisions) is written separately.

    python make_first_sample.py REVIEW_DIR RESOLVED_JSON PACKET_DIR KEY_TSV

The first draw of the sample (its key is first_sample.tsv); make_packets.py extends it to the 400 assertions of key.tsv.
Its inputs are the review's own files (REVIEW_DIR: the test packets as sent to the reviewer, PROMPT.md and
PRO_TEST_REVIEW_*_of_41.md; RESOLVED_JSON: the final decisions), which are not part of the release; the released
benchmark_reviews table holds the same decisions.
"""
import csv
import json
import random
import sys
from pathlib import Path

review_dir, resolved, packet_dir, key = (Path(a) for a in sys.argv[1:5])
SEED, PER_CLASS, PACKETS = 20261030, 100, 4

decisions = {a['assertion_id']: a['decision'] for a in json.load(open(resolved, encoding='utf-8'))}
proteins, assertion = {}, {}
for f in sorted(review_dir.glob('PRO_TEST_REVIEW_*_of_41.md')):
    text = f.read_text(encoding='utf-8')
    batch = json.loads(text.split('```json\n', 1)[1].rsplit('\n```', 1)[0])
    for p in batch['proteins']:
        proteins[p['sequence_sha256']] = {k: v for k, v in p.items() if k != 'assertions'}
        for a in p['assertions']:
            assertion[a['assertion_id']] = (p['sequence_sha256'], a)
if set(assertion) != set(decisions):
    raise ValueError('the packets and the final decisions list different assertions')

rng = random.Random(SEED)
sample = []
for decision in ('accept', 'do_not_accept'):
    pool = sorted(i for i, d in decisions.items() if d == decision)
    sample += [(i, decision) for i in rng.sample(pool, PER_CLASS)]
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
        rows.append(dict(packet=k + 1, assertion_id=i, original_decision=decision))
    data = dict(batch=k + 1, source_release='GOA 2026-07-29',
                proteins=[dict(proteins[seq], assertions=a) for seq, a in grouped.items()])
    name = f'AGREEMENT_RESULTS_{k + 1}.json'
    text = (f'# Protein function review — batch {k + 1} of {PACKETS}\n\n{len(grouped)} proteins; {len(part)} assertions. '
            f'Name your output `{name}`.\n\n{prompt.rstrip()}\n\n\n## Protein and assertion data\n\n```json\n'
            + json.dumps(data, ensure_ascii=False, indent=1) + '\n```\n')
    (packet_dir / f'AGREEMENT_REVIEW_{k + 1}_of_{PACKETS}.md').write_text(text, encoding='utf-8')
    print(f'packet {k + 1}: {len(grouped)} proteins, {len(part)} assertions, '
          f'{sum(d == "accept" for _, d in part)} originally accepted')
with open(key, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader()
    w.writerows(rows)
