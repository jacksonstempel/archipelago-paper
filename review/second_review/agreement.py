"""Agreement between the test-label review (its final decisions, as released in benchmark_reviews) and an independent
second review by a different language model (Claude Opus 5.5, same prompt and packet format, blind to the first;
Supplementary Note 3) on 100 test assertions: packets 1 to 4 of the sample in key.tsv (400 assertions, 200 accepted and 200
not accepted, shuffled into 16 packets of 25 by make_packets.py; the second review was run on the first four packets,
packets/, and returned results/).

Per stratum of the first review: how often the second review agrees. Cohen's kappa on the reviewed sample, and with each
stratum weighted to its share of all 4,429 test assertions (3,538 accepted, 891 not), with 95% intervals from 10,000
resamples of assertions within each stratum. A second-review "needs_adjudication" counts as not accepted (the rubric
accepts only established support) and is reported. Assertions the second review did not finish are left out.

    python agreement.py > agreement.txt

Expected output (Supplementary Note 3): first accept 55, second accepts 51; first do_not_accept 45, second accepts 11;
sample agreement 0.85 [0.78, 0.92], Cohen kappa 0.69 [0.54, 0.84]; population agreement 0.89 [0.83, 0.95], kappa 0.67
[0.51, 0.83].
"""
import csv
import json
from pathlib import Path

import numpy as np

H = Path(__file__).resolve().parent
POPULATION = {'accept': 3538, 'do_not_accept': 891}

key = {r['assertion_id']: r for r in csv.DictReader(open(H / 'key.tsv', encoding='utf-8'), delimiter='\t')}
second, unfinished, adjudication = {}, [], 0
for f in sorted((H / 'results').glob('AGREEMENT_RESULTS_*.json')):
    data = json.load(open(f, encoding='utf-8'))
    for a in data['assertions']:
        if a['assertion_id'] not in key:
            raise ValueError(f'{f.name}: unknown assertion {a["assertion_id"]}')
        adjudication += a['decision'] == 'needs_adjudication'
        second[a['assertion_id']] = 'accept' if a['decision'] == 'accept' else 'do_not_accept'
    unfinished += data.get('unreviewed', [])
print(f'second review: {len(second)} assertions decided ({adjudication} needs_adjudication, counted as not accepted); '
      f'{len(unfinished)} unfinished; {len(key) - len(second) - len(unfinished)} in packets without results')

first = {i: key[i]['original_decision'] for i in second}
both = {s: np.array([second[i] == 'accept' for i in second if first[i] == s]) for s in POPULATION}
for s, x in both.items():
    agree = x.mean() if s == 'accept' else 1 - x.mean()
    print(f'first {s:13s}: {len(x):3d} assertions; second accepts {int(x.sum()):3d}, does not {int((~x).sum()):3d}; '
          f'agreement {agree:.2f}')


def kappa(acc, rej, weights):
    """acc, rej: second-review accept flags in the first review's accepted and not-accepted strata."""
    wa, wr = weights
    a, b = wa * acc.mean(), wa * (1 - acc.mean())  # first accept: second accept / not
    c, d = wr * rej.mean(), wr * (1 - rej.mean())  # first not: second accept / not
    n = a + b + c + d
    po = (a + d) / n
    pe = ((a + b) * (a + c) + (c + d) * (b + d)) / n ** 2
    return po, (po - pe) / (1 - pe)


rng = np.random.default_rng(20261030)
acc, rej = both['accept'], both['do_not_accept']
for name, weights in (('sample', (len(acc), len(rej))), ('population', (POPULATION['accept'], POPULATION['do_not_accept']))):
    po, k = kappa(acc, rej, weights)
    reps = np.array([kappa(rng.choice(acc, len(acc)), rng.choice(rej, len(rej)), weights) for _ in range(10000)])
    lo, hi = np.percentile(reps, [2.5, 97.5], axis=0)  # each: [agreement, kappa]
    print(f'{name:10s}: agreement {po:.2f} [{lo[0]:.2f}, {hi[0]:.2f}], Cohen kappa {k:.2f} [{lo[1]:.2f}, {hi[1]:.2f}]')
