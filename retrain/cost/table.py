"""The arithmetic of the training-cost table (Supplementary Note 9): the published methods' cost with their authors' original
code, from the times the timing harnesses measured (retrain.cost.deepgose, dpfunc, transfew, on one H100 each) and the
times the paper's runs recorded (retrain.scoring.training_cost reads the same records), beside the cost as trained here.

    python -m retrain.cost.table [--out FILE]

The constants below are the values the paper's timing jobs printed and the paper's runs recorded. Hours are H100-hours; one
run is one ensemble of ten worlds for DeepGO-SE and one run (seed) for TransFew and DPFunc; "all" is the four training sets
of the main comparison (original data, HQ with GOA proteins only, HQ, ArchipelaGO), three runs per set for TransFew and
DPFunc. Standard library only.
"""
import argparse
from math import ceil

H = 3600.0

# ---------------------------------------------------------------- DeepGO-SE (retrain.cost.deepgose)
JULY_EPOCH_ARCH = 939.9      # s, 6,427 batches x 145.49 ms + 4.76 s shuffle, the paper's axioms, ArchipelaGO (measured)
MARCH_EPOCH_ARCH = 1009.0    # s, same job, March 2026 axioms (the earlier native run recorded 1,013.0 s: harness check)
RATIO = 145.49 / 156.26      # paper's / March axioms, ms per batch, same job (0.931)
ITEM_MS = 145.85 - 145.49    # official train.py reads both losses with .item() each batch (measured, paper's axioms)
ROC_S = 4.59                 # official compute_roc on the 2,500 development proteins, once per epoch (measured, CPU)
# set: (proteins, epochs run by the paper's ten worlds, March native epoch s measured earlier (GOA: HQ's scaled by
# proteins), hours outside the epochs (recorded))
DEEPGO = {'orig': (35898, 816, 156.4, 0.49), 'goa': (69412, 282, 298.2, 0.19), 'hq': (87797, 261, 377.2, 0.17),
          'archipelago': (237797, 237, None, 0.17)}
DEEPGO_TRAINED = {'orig': 5.03, 'goa': 3.20, 'hq': 3.68, 'archipelago': 8.81}   # the paper's engine, recorded

# ---------------------------------------------------------------- DPFunc (retrain.cost.dpfunc)
AUTH_MS, OURS_MS = 182.9, 212.4     # ms per batch of 64, in-RAM released path vs ours (memory-mapped bank), same GPU, warm
DEV_AUTH, DEV_OURS = 193.95, 7.91   # s per epoch: released test_performance_gnn_inter vs our evaluate, 2,500 dev proteins
DPFUNC = {'orig': (28051, [2661, 2573, 2467]), 'goa': (69412, [8029, 7951, 6026]), 'hq': (87797, [7203, 7369, 7366]),
          'archipelago': (237797, [24489, 12316, 12245])}   # (proteins, recorded training seconds of seeds 0, 1, 2)

# ---------------------------------------------------------------- TransFew (retrain.cost.transfew)
TF_OFF, TF_REL = 18.5, 120.6        # ms per batch of 256 (mean of epochs 1-2): diagnostics off / as released
TF_VAL_OFF, TF_VAL_REL = 8.95, 97.9  # ms per validation batch, 10 batches (2,500 dev proteins)
TF_OUTSIDE = 28.0                   # s per fit outside the epochs (recorded fits)
TRANSFEW = {'orig': (70119, [806, 816, 818]), 'goa': (69412, [782, 786, 809]), 'hq': (87797, [961, 951, 940]),
            'archipelago': (237797, [2265, 2268, 2301])}   # (proteins, recorded GPU seconds of seeds 0, 1, 2)


def table():
    out = []
    P = lambda *a: out.append(' '.join(str(x) for x in a))  # noqa: E731
    SETS = ['orig', 'goa', 'hq', 'archipelago']

    P('DeepGO-SE: one ensemble (ten worlds), native step, the paper\'s epochs and axioms')
    P(f'{"set":12s}{"epochs":>7s}{"epoch s (off)":>14s}{"outside h":>10s}{"off h":>8s}{"released h":>11s}{"trained h":>10s}')
    dg = {}
    for s in SETS:
        n, ep, march, outside = DEEPGO[s]
        epoch = JULY_EPOCH_ARCH if s == 'archipelago' else march * RATIO
        batches = ceil(n / 37)
        rel_epoch = epoch + batches * ITEM_MS / 1000 + ROC_S
        off, rel = (ep * epoch / H + outside), (ep * rel_epoch / H + outside)
        dg[s] = (off, rel)
        P(f'{s:12s}{ep:7d}{epoch:14.1f}{outside:10.2f}{off:8.1f}{rel:11.1f}{DEEPGO_TRAINED[s]:10.2f}')
    P(f'{"all four":12s}{"":7s}{"":14s}{"":10s}{sum(v[0] for v in dg.values()):8.1f}{sum(v[1] for v in dg.values()):11.1f}'
      f'{sum(DEEPGO_TRAINED.values()):10.2f}')
    P('(archipelago epoch = measured with the paper\'s axioms; other sets = March epoch x 0.931; orig and hq measured natively '
      'in March on those sets, goa scaled from hq by proteins)')

    P('\nDPFunc: three runs per set, 15 epochs; recorded time + 15 x (batches x (ours - authors ms) - (authors dev check - ours))')
    P(f'{"set":12s}{"batches":>8s}{"recorded h (3 runs)":>20s}{"authors h (3 runs)":>19s}{"authors h / run":>16s}')
    dp = {}
    for s in SETS:
        n, secs = DPFUNC[s]
        if s == 'archipelago':
            secs = [sum(secs[1:]) / 2] * 3     # seed 0 ran on a slow node (24,489 s); basis: seeds 1 and 2
        b = ceil(n / 64)
        delta = 15 * (b * (AUTH_MS - OURS_MS) / 1000 + (DEV_AUTH - DEV_OURS))
        rec = sum(secs) / H
        dp[s] = (rec, rec + 3 * delta / H, (rec + 3 * delta / H) / 3)
        P(f'{s:12s}{b:8d}{rec:20.2f}{dp[s][1]:19.2f}{dp[s][2]:16.2f}')
    P(f'{"all four":12s}{"":8s}{sum(v[0] for v in dp.values()):20.2f}{sum(v[1] for v in dp.values()):19.2f}')
    P('(archipelago recorded = 3 x mean of seeds 1 and 2; one run as trained here = ', round(dp['archipelago'][0] / 3, 2), 'h)')

    P('\nTransFew: three runs per set, 150 epochs; epoch = batches(256) x per-batch + 10 validation batches; + 28 s per fit')
    P(f'{"set":12s}{"batches":>8s}{"off h/run":>10s}{"released h/run":>15s}{"off h (3)":>10s}{"released h (3)":>15s}{"trained h (3)":>14s}')
    tf = {}
    for s in SETS:
        n, secs = TRANSFEW[s]
        b = ceil(n / 256)
        off = (150 * (b * TF_OFF + 10 * TF_VAL_OFF) / 1000 + TF_OUTSIDE) / H
        rel = (150 * (b * TF_REL + 10 * TF_VAL_REL) / 1000 + TF_OUTSIDE) / H
        tf[s] = (off, rel, sum(secs) / H)
        P(f'{s:12s}{b:8d}{off:10.2f}{rel:15.2f}{3 * off:10.2f}{3 * rel:15.2f}{sum(secs) / H:14.2f}')
    P(f'{"all four":12s}{"":8s}{"":10s}{"":15s}{3 * sum(v[0] for v in tf.values()):10.2f}{3 * sum(v[1] for v in tf.values()):15.2f}'
      f'{sum(v[2] for v in tf.values()):14.2f}')

    P('\nFINAL TABLE (H100-hours): one run on ArchipelaGO / all training sets as in the paper')
    P(f'{"method":10s}{"original, diagnostics off":>28s}{"original, as released":>24s}{"as trained here":>20s}')
    a = 'archipelago'
    P(f'{"DeepGO-SE":10s}{dg[a][0]:>14.1f} / {sum(v[0] for v in dg.values()):<11.1f}{dg[a][1]:>10.1f} / {sum(v[1] for v in dg.values()):<11.1f}'
      f'{DEEPGO_TRAINED[a]:>6.1f} / {sum(DEEPGO_TRAINED.values()):<9.1f}  (one ensemble)')
    P(f'{"DPFunc":10s}{dp[a][2]:>14.1f} / {sum(v[1] for v in dp.values()):<11.1f}{dp[a][2]:>10.1f} / {sum(v[1] for v in dp.values()):<11.1f}'
      f'{dp[a][0] / 3:>6.1f} / {sum(v[0] for v in dp.values()):<9.1f}  (one run / three runs per set)')
    P(f'{"TransFew":10s}{tf[a][0]:>14.2f} / {3 * sum(v[0] for v in tf.values()):<11.1f}{tf[a][1]:>10.2f} / {3 * sum(v[1] for v in tf.values()):<11.1f}'
      f'{tf[a][2] / 3:>6.2f} / {sum(v[2] for v in tf.values()):<9.1f}  (one run / three runs per set)')
    # {method: (original one run, original all, as trained one run, as trained all)}, diagnostics off
    final = {'deepgose': (dg[a][0], sum(v[0] for v in dg.values()), DEEPGO_TRAINED[a], sum(DEEPGO_TRAINED.values())),
             'dpfunc': (dp[a][2], sum(v[1] for v in dp.values()), dp[a][0] / 3, sum(v[0] for v in dp.values())),
             'transfew': (tf[a][0], 3 * sum(v[0] for v in tf.values()), tf[a][2] / 3, sum(v[2] for v in tf.values()))}
    return '\n'.join(out) + '\n', final


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--out', help='also write the table to this file')
    a = parser.parse_args()
    text, _ = table()
    if a.out:
        with open(a.out, 'w', encoding='utf-8', newline='\n') as handle:
            handle.write(text)
    print(text, end='')


if __name__ == '__main__':
    main()
