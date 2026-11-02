"""Figure 2 (every published method retrained), Figure 4 (Finch: the recipe ablation, and Finch against the published
methods on every training set they share) and the bodies of Tables 1 and 2, from the results: <root>/results/main/
{metrics,contrasts}.tsv (test, seed means), <root>/results/explabels_training/{metrics,contrasts}.tsv and <root>/results/
ladder/metrics.tsv; Table 2 also reads the contrasts.tsv of every version of the test set (main,
experimental_test_labels, masked_test_labels, all_goa_test_labels, no_hit/E0.001, no_hit/E10) with the label versions'
counts.tsv and no_hit/subsets.tsv."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import blended_transform_factory

from figures.style import COLOURS, DOUBLE, GREY, NAMES, PUBLISHED, contrasts, fs, half_up, read, save, seed_means

FINCH = COLOURS['archipelago']
# Fig. 4a: the recipe ablation's steps (configs/finch.json, 'ladder'), ending at Finch; its control (two ESM2-15B networks)
# is not drawn
LADDER = [('ladder_linear_esm2_3b', 'Linear readout,\nESM2-3B'), ('ladder_linear_esm2_15b', 'ESM2-15B'),
          ('ladder_narrow', '+ hidden\nlayer'), ('ladder_hidden', '+ width,\nheavy\ndropout'), ('ladder_asam', '+ ASAM'),
          ('finch', '+ ProtT5\nnetwork')]
Y0, Y1 = 0.45, 0.65  # Fig. 4a's axis
HALO = dict(facecolor='white', edgecolor='none', pad=0.5)  # behind text a method's dashed line passes through
MARK = {'deepgose': ('o', 4.2), 'transfew': ('s', 3.8), 'dpfunc': ('^', 4.6)}  # Fig. 4b: marker and size of each method
DOT_INK = '#55595f'  # Fig. 4b: the published methods


def bar(ax, x, v, color, w):
    ax.bar(x, v, w, color=color, lw=0, zorder=2)
    ax.text(x, v + 0.01, half_up(v), ha='center', va='bottom', fontsize=fs(6.5), zorder=4)


def style_bars(ax):
    ax.set_ylim(0, 0.7)
    ax.set_yticks(np.arange(0, 0.71, 0.1))
    ax.set_ylabel('Micro-AUPRC')
    ax.tick_params(axis='x', length=2.5, width=0.6)


def waterfall(ax, auprc):
    """Fig. 4a: test micro-AUPRC from a linear readout on ESM2-3B (a full bar) through one ingredient per step (floating
    bars from the level before) to Finch (a full bar), each label under its bar; the published methods as dashed lines. The
    hidden layer enters in two steps: narrow (1,024 units) with light dropout, then wide (8,192) with heavy dropout."""
    levels = [auprc(m, 'archipelago') for m, _ in LADDER]
    w, gap = 0.62, 0.004  # bar width; room between a bar's bottom and its label
    end = len(levels)  # Finch, the full bar after the last step
    below_axis = blended_transform_factory(ax.transData, ax.transAxes)

    def label(x, y, text):
        """A step's name under its bar: under the axis for the full bars, under the bar's bottom for the steps."""
        if y is None:
            ax.text(x, -0.02, text, transform=below_axis, ha='center', va='top', fontsize=fs(6.5), linespacing=1.05)
        else:
            ax.text(x, y - gap, text, ha='center', va='top', fontsize=fs(6.5), linespacing=1.05, zorder=4, bbox=HALO)

    ax.bar(0, levels[0] - Y0, w, bottom=Y0, color='#8f939a', lw=0, zorder=3)  # the start, a full bar
    ax.text(0, levels[0] + 0.003, half_up(levels[0]), ha='center', va='bottom', fontsize=fs(6.5), zorder=4, bbox=HALO)
    label(0, None, LADDER[0][1])
    for i in range(1, len(levels)):
        lo, hi = levels[i - 1], levels[i]
        ax.bar(i, hi - lo, w, bottom=lo, color=FINCH, lw=0, zorder=3)
        ax.plot([i - 1 + w / 2, i - w / 2], [lo, lo], color=GREY, lw=0.6, ls=(0, (1.5, 1.5)), zorder=2)  # connector
        ax.text(i, hi + 0.003, half_up(hi - lo, sign=True), ha='center', va='bottom', fontsize=fs(6.5), zorder=4, bbox=HALO)
        label(i, lo, LADDER[i][1])
    ax.plot([end - 1 + w / 2, end - w / 2], [levels[-1]] * 2, color=GREY, lw=0.6, ls=(0, (1.5, 1.5)), zorder=2)
    ax.bar(end, levels[-1] - Y0, w, bottom=Y0, color=FINCH, lw=0, zorder=3)  # Finch, a full bar
    ax.text(end, levels[-1] + 0.003, half_up(levels[-1]), ha='center', va='bottom', fontsize=fs(6.5), zorder=4)
    label(end, None, 'Finch')
    # the published methods trained on ArchipelaGO, dashed behind the bars. Their labels stand in a column against the Finch
    # bar, the one stretch free of step labels; DeepGO-SE and TransFew are 0.001 apart, so one sits above their lines and the
    # other below
    for m, va, dy in (('deepgose', 'bottom', 0.0012), ('transfew', 'top', -0.0012), ('dpfunc', 'bottom', 0.0012)):
        v = auprc(m, 'archipelago')
        ax.axhline(v, color='#8f939a', lw=0.8, ls=(0, (3, 2)), zorder=1)
        ax.text(end - w / 2 - 0.08, v + dy, f'{NAMES[m]} {half_up(v)}', ha='right', va=va, fontsize=fs(6.5), color=GREY, zorder=4)
    ax.spines['bottom'].set_zorder(5)  # the axis line over the bars' bottoms
    ax.set_xticks([])
    ax.set_xlim(-0.55, end + 0.55)
    ax.set_ylim(Y0, Y1)
    ax.set_yticks(np.arange(Y0, Y1 + 0.001, 0.05))
    ax.set_ylabel('Micro-AUPRC')


def shared_rows(results):
    """Fig. 4b's rows, lowest Finch score first: (label, {published method: micro-AUPRC}, Finch's, Finch's lead over the
    row's best published method), each published method's original data first."""
    main, trained = seed_means(results / 'main'), seed_means(results / 'explabels_training')
    main_pairs, trained_pairs = contrasts(results / 'main'), contrasts(results / 'explabels_training')
    rows = []
    for m in PUBLISHED:
        own = f'original_{m}'
        rows.append((f"{NAMES[m]}'s\noriginal data", {m: float(main[m, 'original', 'auprc']['value'])},
                     float(main['finch', own, 'auprc']['value']),
                     float(main_pairs['auprc', 'finch', own, m, 'original']['difference'])))
    for s, label, metric, pairs in (('hq_explabels', 'HQ, experimental\nlabels only', trained, trained_pairs),
                                    ('goa', 'HQ, GOA\nproteins only', main, main_pairs), ('hq', 'HQ', main, main_pairs),
                                    ('archipelago', 'ArchipelaGO', main, main_pairs)):
        published = {m: float(metric[m, s, 'auprc']['value']) for m in PUBLISHED}
        best = max(published, key=published.get)
        rows.append((label, published, float(metric['finch', s, 'auprc']['value']),
                     float(pairs['auprc', 'finch', s, best, s]['difference'])))
    return sorted(rows, key=lambda r: r[2])


def dots(ax, rows, comparators):
    """Fig. 4b: one row per training set, the published methods' dots, Finch's dot with its score, a line from the row's
    best published method to Finch, and Finch's lead on the right."""
    n = len(rows)
    for i, (_, published, finch, gain) in enumerate(rows):
        y = n - 1 - i
        ax.axhline(y, color='#e6e7e9', lw=0.6, zorder=1)
        ax.plot([max(published.values()), finch], [y, y], color=FINCH, lw=1.4, alpha=0.35, solid_capstyle='butt', zorder=2)
        for m, v in published.items():
            near = [o for o in published if o != m and abs(published[o] - v) < 0.004]  # would sit on top of each other
            shift = 0 if not near else (0.17 if v > published[near[0]] else -0.17)
            ax.plot(v, y + shift, MARK[m][0], ms=MARK[m][1], mfc=DOT_INK, mec='white', mew=0.5, zorder=3)
        ax.plot(finch, y, 'o', ms=5.6, mfc=FINCH, mec='white', mew=0.5, zorder=4)
        ax.text(finch, y + 0.2, half_up(finch), ha='center', va='bottom', fontsize=fs(6.5), zorder=5)
        ax.text(1.02, y, half_up(gain, sign=True), transform=ax.get_yaxis_transform(), ha='left', va='center', fontsize=fs(6.5))
    ax.set_yticks(range(n), [r[0] for r in rows][::-1], fontsize=fs(6.5), linespacing=1.05)
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    ax.set_ylim(-0.6, n - 0.35)
    ax.set_xlim(0.38, 0.64)
    ax.set_xticks(np.arange(0.40, 0.601, 0.05))
    ax.set_xlabel('Micro-AUPRC')
    handles = [plt.Line2D([], [], ls='', marker='o', ms=5.6, mfc=FINCH, mec='white', mew=0.5, label='Finch')]
    handles += [plt.Line2D([], [], ls='', marker=MARK[m][0], ms=MARK[m][1], mfc=DOT_INK, mec='white', mew=0.5, label=NAMES[m])
                for m in comparators]
    ax.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=4, frameon=False, fontsize=fs(6.5),
              handletextpad=0.0, columnspacing=0.6, borderaxespad=0.1, borderpad=0)


def figures(results, out):
    metric = seed_means(results / 'main')
    auprc = lambda m, s: float(metric[m, s, 'auprc']['value'])  # noqa: E731
    comparators = sorted(PUBLISHED, key=lambda m: -auprc(m, 'original'))  # Fig. 2 and Table 1 order

    # Figure 4. a: the recipe ablation as a waterfall; b: Finch and the published methods on every training set they share.
    ladder = seed_means(results / 'ladder')
    ladder_auprc = lambda m, s: float(ladder[m, s, 'auprc']['value'])  # noqa: E731
    fig = plt.figure(figsize=(DOUBLE, 3.15))
    grid = fig.add_gridspec(1, 2, width_ratios=[1.3, 1], wspace=0.42, left=0.06, right=0.945, top=0.91, bottom=0.14)
    a, b = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
    waterfall(a, ladder_auprc)
    dots(b, shared_rows(results), comparators)
    for ax, letter, x in ((a, 'a', -0.11), (b, 'b', -0.3)):
        ax.text(x, 1.02, letter, transform=ax.transAxes, fontsize=8, fontweight='bold', va='bottom', ha='left')
    save(fig, out, 'fig4_finch')

    # Figure 2. Per published method, original data, HQ with experimental labels only, HQ and ArchipelaGO side by side;
    # methods by original-data score.
    W = 0.25  # one bar
    GAP = 0.25 * 0.85  # between the four-bar groups
    order = [('original', COLOURS['original'], 'Original data'), ('hq_explabels', '#aad2ef', 'HQ, experimental labels only'),
             ('hq', COLOURS['hq'], 'HQ'), ('archipelago', FINCH, 'ArchipelaGO')]
    scores = {**metric, **{k: r for k, r in seed_means(results / 'explabels_training').items() if k[1] == 'hq_explabels'}}
    bars = len(order)
    fig, ax = plt.subplots(figsize=(DOUBLE * 0.9, 2.55))
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.11, top=0.88)  # fixed, so the key can be placed against the axes
    centres = []
    x = bars / 2 * W  # half a group from the axis
    for m in comparators:
        centres.append(x)
        for k, (s, color, _) in enumerate(order):
            bar(ax, x + (k - (bars - 1) / 2) * W, float(scores[m, s, 'auprc']['value']), color, W * 0.94)
        x += bars * W + GAP
    ax.set_xticks(centres, [NAMES[m] for m in comparators], fontsize=fs(7))
    ax.set_xlim(-0.15, centres[-1] + bars / 2 * W + 0.15)
    style_bars(ax)
    # the key: 'Trained on:' above the first group, the swatches right after it
    lead = ax.text(0.0, 1.06, 'Trained on:', transform=blended_transform_factory(ax.transData, ax.transAxes), ha='left',
                   va='center', fontsize=fs(6.5))
    fig.canvas.draw()
    right = ax.transAxes.inverted().transform(lead.get_window_extent().get_points())[1, 0]
    ax.legend([plt.Rectangle((0, 0), 1, 1, color=c) for _, c, _ in order], [t for *_, t in order], loc='center left',
              bbox_to_anchor=(right + 0.008, 1.06), ncol=bars, frameon=False, fontsize=fs(6.5), handlelength=1.8,
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


def table2(results, out):
    """Table 2: the paper's three comparisons in test micro-AUPRC (each published method's gain from its original data to
    ArchipelaGO; Finch on ArchipelaGO minus Finch on HQ + random, and minus each published method on ArchipelaGO), on
    every version of the test set, read from each version's contrasts.tsv. table2.tex is the tabular the manuscript
    inputs (full text width, the reviewed labels shaded: colortbl, with `rowshade` defined by the document);
    table2_intervals.tex, the SI's copy at natural width, adds the paired 95% cluster-bootstrap intervals."""
    metric = seed_means(results / 'main')
    comparators = sorted(PUBLISHED, key=lambda m: -float(metric[m, 'original', 'auprc']['value']))
    proteins = lambda folder, variant: next(int(r['proteins']) for r in read(results / folder / 'counts.tsv')  # noqa: E731
                                            if r['variant'] == variant)
    subsets = {float(r['cutoff']): int(r['proteins']) for r in read(results / 'no_hit' / 'subsets.tsv')}
    # (name, proteins, results folder, main result) per group: the label versions from strictest to loosest, then the
    # protein subsets; the groups are separated by a rule, not labelled
    groups = [[('Experimental labels only', proteins('experimental_test_labels', 'exp'), 'experimental_test_labels', False),
               ('Annotations masked', proteins('masked_test_labels', 'masked'), 'masked_test_labels', False),
               ('Reviewed labels', proteins('experimental_test_labels', 'all'), 'main', True),
               ('All GOA annotations', proteins('all_goa_test_labels', 'all_goa'), 'all_goa_test_labels', False)],
              [(r'No training hit, E\,$\leq$\,10\textsuperscript{-3}', subsets[1e-3], 'no_hit/E0.001', False),
               (r'No training hit, E\,$\leq$\,10', subsets[10.0], 'no_hit/E10', False)]]
    n = len(comparators)
    for name, intervals in (('table2.tex', False), ('table2_intervals.tex', True)):
        # full width through tabularx's stretched first column, not tabular*'s stretched gaps, which a shaded row leaves white
        lines = [r'\begin{tabular}{@{}lr' + 'c' * (2 * n + 1) + '@{}}' if intervals else
                 r'\begin{tabularx}{\textwidth}{@{}Xr' + 'c' * (2 * n + 1) + '@{}}',
                 r'\toprule',
                 rf' & & \multicolumn{{{n}}}{{c}}{{\bfseries Original data $\rightarrow$ ArchipelaGO}} & '
                 rf'\multicolumn{{{n + 1}}}{{c}}{{\bfseries Finch on ArchipelaGO minus}} \\',
                 rf'\cmidrule(lr){{3-{2 + n}}}\cmidrule(lr){{{3 + n}-{3 + 2 * n}}}',
                 r'\bfseries Test set & Proteins & ' + ' & '.join(NAMES[m] for m in comparators) + r' & HQ + random & '
                 + ' & '.join(NAMES[m] for m in comparators) + r' \\']
        for group in groups:
            lines.append(r'\midrule')
            for version, count, folder, main in group:
                c = contrasts(results / folder)
                row = [c['auprc', m, 'archipelago', m, 'original'] for m in comparators] + \
                      [c['auprc', 'finch', 'archipelago', 'finch', 'hq_random']] + \
                      [c['auprc', 'finch', 'archipelago', m, 'archipelago'] for m in comparators]
                cells = [f"{half_up(r['difference'])} ({half_up(r['ci_low'])}--{half_up(r['ci_high'])})" if intervals
                         else half_up(r['difference']) for r in row]
                lines.append((r'\rowcolor{rowshade}[0pt][0pt]' if main else '') + f'{version} & {count:,} & ' + ' & '.join(cells)
                             + r' \\')
        lines += [r'\bottomrule', r'\end{tabular}' if intervals else r'\end{tabularx}', '']
        with open(out / name, 'w', encoding='utf-8', newline='\n') as handle:
            handle.write('\n'.join(lines))
