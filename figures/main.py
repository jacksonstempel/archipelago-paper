"""Figure 2 (every published method retrained), Figure 4 (Finch: the recipe ablation and precision-recall curves) and the
body of Table 1, from the results: <root>/results/main/metrics.tsv (test, seed means), <root>/results/ladder/metrics.tsv
and <root>/results/pr_curves/curves.tsv."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import blended_transform_factory

from figures.style import COLOURS, DOUBLE, GREY, INK, NAMES, PUBLISHED, fs, half_up, pr_curves, save, seed_means

FINCH = COLOURS['archipelago']
SHADES = ['#55595f', '#8f939a', '#bfc2c7']  # the published methods, best on ArchipelaGO first
# Fig. 4a: the recipe ablation's steps (configs/finch.json, 'ladder'), ending at Finch; its control (two ESM2-15B networks)
# is not drawn
LADDER = [('ladder_linear_esm2_3b', 'Linear readout\nESM2-3B'), ('ladder_linear_esm2_15b', 'ESM2-15B'),
          ('ladder_hidden', '+ Hidden\nlayer'), ('ladder_asam', '+ ASAM'), ('finch', '+ ProtT5\nnetwork (Finch)')]
SHARE = 0.002  # published methods closer than this to the one above them share its label in Fig. 4a


def bar(ax, x, v, color, w):
    ax.bar(x, v, w, color=color, lw=0, zorder=2)
    ax.text(x, v + 0.01, half_up(v), ha='center', va='bottom', fontsize=fs(6.5), zorder=4)


def style_bars(ax):
    ax.set_ylim(0, 0.7)
    ax.set_yticks(np.arange(0, 0.71, 0.1))
    ax.set_ylabel('Micro-AUPRC')
    ax.tick_params(axis='x', length=2.5, width=0.6)


def waterfall(ax, auprc):
    """Fig. 4a: test micro-AUPRC from a linear readout on ESM2-3B to Finch, one ingredient per step, each step a floating bar
    from the level before it; the published methods trained on ArchipelaGO as dashed lines."""
    levels = [auprc(m, 'archipelago') for m, _ in LADDER]
    w = 0.62
    # the start: a level, not a bar (the axis does not start at zero)
    ax.plot([-w / 2, w / 2], [levels[0]] * 2, color=GREY, lw=2.2, solid_capstyle='butt', zorder=3)
    ax.text(0, levels[0] - 0.003, half_up(levels[0]), ha='center', va='top', fontsize=fs(6.5), zorder=4)
    for i in range(1, len(levels)):
        lo, hi = levels[i - 1], levels[i]
        ax.bar(i, hi - lo, w, bottom=lo, color=FINCH, lw=0, zorder=3)
        ax.plot([i - 1 + w / 2, i - w / 2], [lo, lo], color=GREY, lw=0.6, ls=(0, (1.5, 1.5)), zorder=2)  # connector
        ax.text(i, hi + 0.0015, half_up(hi - lo, sign=True), ha='center', va='bottom', fontsize=fs(6.5), zorder=4)
    ax.text(len(levels) - 1 + w / 2 + 0.06, levels[-1], half_up(levels[-1]), ha='left', va='center', fontsize=fs(6.5),
            fontweight='bold', zorder=4)
    labels = []  # [value, names]: a method close to the one above it shares that one's label
    for m in sorted(PUBLISHED, key=lambda m: -auprc(m, 'archipelago')):
        v = auprc(m, 'archipelago')
        ax.axhline(v, color='#8f939a', lw=0.8, ls=(0, (3, 2)), zorder=1)
        if labels and labels[-1][0] - v < SHARE:
            labels[-1][1].append(NAMES[m])
        else:
            labels.append([v, [NAMES[m]]])
    for v, names in labels:
        ax.text(-0.45, v + 0.0012, ', '.join(names), ha='left', va='bottom', fontsize=fs(6.5), color=GREY)
    ax.set_xticks(range(len(LADDER)), [t for _, t in LADDER], fontsize=fs(6.5), linespacing=1.05)
    ax.tick_params(axis='x', length=0)
    ax.set_xlim(-0.55, len(LADDER) - 0.3)
    ax.set_ylim(0.52, 0.64)
    ax.set_yticks(np.arange(0.52, 0.641, 0.02))
    ax.set_ylabel('Micro-AUPRC')


def figures(results, out):
    metric = seed_means(results / 'main')
    auprc = lambda m, s: float(metric[m, s, 'auprc']['value'])  # noqa: E731
    comparators = sorted(PUBLISHED, key=lambda m: -auprc(m, 'original'))  # Fig. 2 and Table 1 order

    # Figure 4. a: the recipe ablation as a waterfall; b: Finch and the published methods, all trained on ArchipelaGO, as
    # micro-averaged precision-recall curves (averaged over runs at each recall, so each area is the model's micro-AUPRC),
    # each with its Fmax point.
    ladder = seed_means(results / 'ladder')
    ladder_auprc = lambda m, s: float(ladder[m, s, 'auprc']['value'])  # noqa: E731
    fig = plt.figure(figsize=(DOUBLE, 3.15))
    grid = fig.add_gridspec(1, 2, width_ratios=[1.3, 1], wspace=0.2, left=0.06, right=0.985, top=0.93, bottom=0.14)
    a, b = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
    waterfall(a, ladder_auprc)
    pr = pr_curves(results)
    handles = []
    lines = [(('finch', 'archipelago'), FINCH, 1.5, 'Finch')] + \
            [((m, 'archipelago'), shade, 1.1, NAMES[m])
             for m, shade in zip(sorted(PUBLISHED, key=lambda m: -auprc(m, 'archipelago')), SHADES)]
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
    for ax, letter, x in ((a, 'a', -0.11), (b, 'b', -0.14)):
        ax.text(x, 1.02, letter, transform=ax.transAxes, fontsize=fs(8), fontweight='bold', va='bottom', ha='left')
    save(fig, out, 'fig4_finch')

    # Figure 2. Per published method, original data, HQ and ArchipelaGO side by side; methods by original-data score.
    W = 0.25  # one bar
    GAP = 0.25 * 0.85  # between the three-bar groups
    order = [('original', COLOURS['original'], 'Original data'), ('hq', COLOURS['hq'], 'HQ'),
             ('archipelago', FINCH, 'ArchipelaGO')]
    fig, ax = plt.subplots(figsize=(DOUBLE * 0.9, 2.8))
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.11, top=0.88)  # fixed, so the key can be placed against the axes
    centres = []
    x = 1.5 * W  # half a group from the axis
    for m in comparators:
        centres.append(x)
        for k, (s, color, _) in enumerate(order):
            bar(ax, x + (k - 1) * W, auprc(m, s), color, W * 0.94)
        x += 3 * W + GAP
    ax.set_xticks(centres, [NAMES[m] for m in comparators], fontsize=fs(7))
    ax.set_xlim(-0.15, centres[-1] + 1.5 * W + 0.15)
    style_bars(ax)
    # the key: 'Trained on:' above the first group, the swatches right after it
    lead = ax.text(0.0, 1.06, 'Trained on:', transform=blended_transform_factory(ax.transData, ax.transAxes), ha='left',
                   va='center', fontsize=fs(6.5))
    fig.canvas.draw()
    right = ax.transAxes.inverted().transform(lead.get_window_extent().get_points())[1, 0]
    ax.legend([plt.Rectangle((0, 0), 1, 1, color=c) for _, c, _ in order], [t for *_, t in order], loc='center left',
              bbox_to_anchor=(right + 0.008, 1.06), ncol=3, frameon=False, fontsize=fs(6.5), handlelength=1.8,
              handleheight=0.45, handletextpad=0.4, columnspacing=1.3, borderaxespad=0, borderpad=0)
    save(fig, out, 'fig2_retraining')

    # Table 1 (the tabular that the manuscript inputs): per model, micro-AUPRC, Fmax and Smin on three training sets.
    sets = [('original', 'Original data'), ('hq', 'HQ'), ('archipelago', 'ArchipelaGO')]
    groups = [['finch'] + comparators, ['embedding_knn', 'diamond_knn', 'naive']]

    def text(m, s, k):
        r = metric.get((m, s, k))
        return None if r is None else (f"{float(r['value']):.2f}" if k == 'smin' else half_up(r['value']))

    # the best value of each metric in the whole table, as printed, is set in bold (lowest Smin, highest AUPRC and Fmax)
    shown = {k: [float(t) for m in sum(groups, []) for s, _ in sets if (t := text(m, s, k)) is not None]
             for k in ('auprc', 'fmax', 'smin')}
    best = {k: min(v) if k == 'smin' else max(v) for k, v in shown.items()}

    def cell(m, s, k):
        t = text(m, s, k)
        return '--' if t is None else (rf'\textbf{{{t}}}' if float(t) == best[k] else t)

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
