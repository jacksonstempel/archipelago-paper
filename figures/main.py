"""Figure 1 (Finch beside the published methods), Figure 3 (every published method retrained) and the body of Table 1, from
the main run's results: <root>/results/main/metrics.tsv (test, seed means) and <root>/results/pr_curves/curves.tsv."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import blended_transform_factory

from figures.style import COLOURS, DOUBLE, INK, NAMES, PUBLISHED, fs, half_up, pr_curves, save, seed_means

FINCH = COLOURS['archipelago']
SHADES = ['#55595f', '#8f939a', '#bfc2c7']  # the published methods, in the order of Fig. 1a


def bar(ax, x, v, color, w):
    ax.bar(x, v, w, color=color, lw=0, zorder=2)
    ax.text(x, v + 0.01, half_up(v), ha='center', va='bottom', fontsize=fs(6.5), zorder=4)


def style_bars(ax):
    ax.set_ylim(0, 0.7)
    ax.set_yticks(np.arange(0, 0.71, 0.1))
    ax.set_ylabel('Micro-AUPRC')
    ax.tick_params(axis='x', length=2.5, width=0.6)


def figures(results, out):
    metric = seed_means(results / 'main')
    auprc = lambda m, s: float(metric[m, s, 'auprc']['value'])  # noqa: E731
    comparators = sorted(PUBLISHED, key=lambda m: -auprc(m, 'original'))

    # Figure 1. a: Finch (ArchipelaGO) beside the published methods on their original data, ordered by score; b: the same
    # comparison as micro-averaged precision-recall curves (averaged over runs at each recall, so each area is the bar's
    # micro-AUPRC), each with its Fmax point.
    fig = plt.figure(figsize=(DOUBLE, 3.15))
    grid = fig.add_gridspec(1, 2, width_ratios=[1.05, 1], wspace=0.2, left=0.075, right=0.985, top=0.93, bottom=0.12)
    a, b = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
    bar(a, 0, auprc('finch', 'archipelago'), FINCH, 0.68)
    for i, m in enumerate(comparators, 1):
        bar(a, i, auprc(m, 'original'), COLOURS['original'], 0.68)
    a.set_xticks(range(4), ['Finch'] + [NAMES[m] for m in comparators], fontsize=fs(7))
    a.set_xlim(-0.6, 3.6)
    style_bars(a)
    pr = pr_curves(results)
    handles = []
    lines = [(('finch', 'archipelago'), FINCH, 1.5, 'Finch')] + \
            [((m, 'original'), shade, 1.1, NAMES[m]) for m, shade in zip(comparators, SHADES)]
    for (m, s), color, lw, name in lines:
        recall, precision = np.array(pr[m, s]).T
        fmax = half_up(metric[m, s, 'fmax']['value'])
        handles.append(b.plot(recall, precision, color=color, lw=lw, solid_capstyle='butt', zorder=3 if m == 'finch' else 2,
                              label=f'{name}  {half_up(auprc(m, s))} | {fmax}')[0])
        f1 = 2 * recall * precision / np.maximum(recall + precision, 1e-12)
        k = int(np.argmax(f1))  # the Fmax point on the curve
        b.scatter(recall[k], precision[k], s=16, color=color, edgecolors=INK, linewidths=0.5, zorder=4)
    b.legend(handles=handles, loc='lower left', title='AUPRC | Fmax', frameon=False, fontsize=fs(6.5),
             title_fontsize=fs(6.5), handlelength=1.6, alignment='left', borderaxespad=0.2)
    b.set_xlim(0, 1)
    b.set_ylim(0, 1.02)
    b.set_xticks(np.arange(0, 1.01, 0.2))
    b.set_yticks(np.arange(0, 1.01, 0.2))
    b.set_box_aspect(1)
    b.set_xlabel('Recall')
    b.set_ylabel('Precision')
    for ax, letter in ((a, 'a'), (b, 'b')):
        ax.text(-0.14, 1.02, letter, transform=ax.transAxes, fontsize=fs(8), fontweight='bold', va='bottom', ha='left')
    save(fig, out, 'fig1_finch')

    # Figure 3. Per published method, ArchipelaGO, HQ and original data side by side, beside Finch; order as Fig. 1.
    W = 0.25  # one bar
    GAP = 0.25 * 0.85  # between the three-bar groups
    order = [('archipelago', FINCH, 'ArchipelaGO'), ('hq', COLOURS['hq'], 'HQ'),
             ('original', COLOURS['original'], 'Original data')]
    fig, ax = plt.subplots(figsize=(DOUBLE * 0.9, 2.8))
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.11, top=0.88)  # fixed, so the key can be placed against the axes
    bar(ax, 0, auprc('finch', 'archipelago'), FINCH, W * 0.94)
    centres = [0.0]
    x = W / 2 + 0.25 * 0.75 + 1.5 * W  # Finch's edge, the gap to the first group (a quarter narrower), half a group
    for m in comparators:
        centres.append(x)
        for k, (s, color, _) in enumerate(order):
            bar(ax, x + (k - 1) * W, auprc(m, s), color, W * 0.94)
        x += 3 * W + GAP
    ax.set_xticks(centres, ['Finch'] + [NAMES[m] for m in comparators], fontsize=fs(7))
    ax.set_xlim(-0.35 * 0.8, centres[-1] + 1.5 * W + 0.15)
    style_bars(ax)
    # the key: 'Trained on:' starting halfway between the y axis and Finch's bar, the swatches right after it
    x0 = (ax.get_xlim()[0] - W * 0.94 / 2) / 2
    lead = ax.text(x0, 1.06, 'Trained on:', transform=blended_transform_factory(ax.transData, ax.transAxes), ha='left',
                   va='center', fontsize=fs(6.5))
    fig.canvas.draw()
    right = ax.transAxes.inverted().transform(lead.get_window_extent().get_points())[1, 0]
    ax.legend([plt.Rectangle((0, 0), 1, 1, color=c) for _, c, _ in order], [t for *_, t in order], loc='center left',
              bbox_to_anchor=(right + 0.008, 1.06), ncol=3, frameon=False, fontsize=fs(6.5), handlelength=1.8,
              handleheight=0.45, handletextpad=0.4, columnspacing=1.3, borderaxespad=0, borderpad=0)
    save(fig, out, 'fig3_retraining')

    # Table 1 (the tabular that the manuscript inputs): per model, micro-AUPRC, Fmax and Smin on three training sets.
    sets = [('original', 'Original data'), ('hq', 'HQ'), ('archipelago', 'ArchipelaGO')]
    groups = [['finch'] + comparators, ['embedding_knn', 'diamond_knn', 'naive']]

    def cell(m, s, k):
        r = metric.get((m, s, k))
        return '--' if r is None else (f"{float(r['value']):.2f}" if k == 'smin' else half_up(r['value']))

    lines = [r'\begin{tabular*}{\textwidth}{@{}l@{\extracolsep{\fill}}' + 'ccc' * len(sets) + '@{}}', r'\toprule',
             ' & ' + ' & '.join(rf'\multicolumn{{3}}{{c}}{{\bfseries {t}}}' for _, t in sets) + r' \\',
             ''.join(rf'\cmidrule(lr){{{2 + 3 * i}-{4 + 3 * i}}}' for i in range(len(sets))),
             r'\bfseries Model' + r' & AUPRC & Fmax & Smin' * len(sets) + r' \\', r'\midrule']
    for g, group in enumerate(groups):
        if g:
            lines.append(r'\midrule')
        for m in group:
            lines.append(f'{NAMES[m]} & ' + ' & '.join(cell(m, s, k) for s, _ in sets for k in ('auprc', 'fmax', 'smin'))
                         + r' \\')
    lines += [r'\bottomrule', r'\end{tabular*}', '']
    with open(out / 'table1.tex', 'w', encoding='utf-8', newline='\n') as handle:
        handle.write('\n'.join(lines))
