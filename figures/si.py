"""Tables and figures of the Supplementary Information: LaTeX table bodies (rows only; the SI typesets the headers) and
two figures, written into <out>/si/. Reads, under <root>/results/:

  main/{metrics.tsv,contrasts.tsv}                         every model on every training set (retrain.scoring.score main)
  dataset_sizes.tsv                                        retrain.scoring.sizes
  experimental_test_labels/counts.tsv                      retrain.scoring.test_labels experimental
  explabels_training/{metrics.tsv,contrasts.tsv}           retrain.scoring.score explabels_training
  no_hit/{subsets.tsv,differences.tsv,E<cutoff>/}          retrain.scoring.no_hit
  pr_curves/curves.tsv                                     retrain.scoring.pr_curves
  selection/contrasts.tsv                                  retrain.scoring.score selection
  ladder/{metrics.tsv,contrasts.tsv}                       retrain.scoring.score ladder
  parameters.json                                          retrain.scoring.parameters
  training_cost.tsv                                        retrain.scoring.training_cost

Main-text Table 2 with intervals (table2_intervals.tex) is written beside Table 2 by figures/main.py.

Supplementary Table 1 (the annotation sources) is typed in the SI and the label-review prompt is released with the
reviewed labels; neither is generated here.
"""
import json

import matplotlib.pyplot as plt
import numpy as np

from figures.style import COLOURS, DOUBLE, INK, NAMES, PUBLISHED, contrasts, fs, half_up, pr_curves, read, save, seed_means
from retrain.finch import config as finch_config

SETS = [('original', 'Original'), ('hq_explabels', 'HQ, experimental labels only'), ('goa', 'HQ, GOA proteins only'), ('hq', 'HQ'),
        ('hq_random', 'HQ + random'), ('archipelago', 'ArchipelaGO')]
SETNAME = dict(SETS)
DASH = '--'
# every training set Finch shares with the published methods, in the order of main-text Fig. 4b: (set, name, the method
# whose original data it is)
SHARED = [('original_dpfunc', "DPFunc's original data", 'dpfunc'), ('original_transfew', "TransFew's original data", 'transfew'),
          ('original_deepgose', "DeepGO-SE's original data", 'deepgose'), ('hq_explabels', 'HQ, experimental labels only', None),
          ('goa', 'HQ, GOA proteins only', None), ('hq', 'HQ', None), ('archipelago', 'ArchipelaGO', None)]
SHADES = {'transfew': '#55595f', 'deepgose': '#8f939a', 'dpfunc': '#bfc2c7'}
# Supplementary Figure 2: an original set wears its method's shade, the HQ family tints of the HQ blue
SET_LINE = {'original_dpfunc': SHADES['dpfunc'], 'original_transfew': SHADES['transfew'], 'original_deepgose': SHADES['deepgose'],
            'hq_explabels': '#aad2ef', 'goa': '#75b6e4', 'hq': COLOURS['hq'], 'archipelago': COLOURS['archipelago']}


def r3(value, sign=False, places=3):
    """Rounded half up, as in the main text, with a LaTeX minus sign."""
    return half_up(value, sign, places).replace('-', '$-$')


def ci(low, high, places=3):
    return f'[{r3(low, places=places)}, {r3(high, places=places)}]'


def bound(value, places):
    """Rounded half up; a negative number that rounds to zero keeps its sign, so an interval shows on which side of zero
    it ends."""
    s = r3(value, places=places)
    return '$-$' + s if float(value) < 0 and not s.startswith('$-$') else s


def paired(d, places):
    """A contrasts.tsv row: the signed difference and its interval."""
    s = bound(d['difference'], places)
    return f"{s if s.startswith('$-$') else '+' + s} [{bound(d['ci_low'], places)}, {bound(d['ci_high'], places)}]"


def write(out, name, lines):
    with open(out / name, 'w', encoding='utf-8', newline='\n') as handle:
        handle.write('\n'.join(lines) + '\n')


def tables(results, out):
    """Every table body of the SI but those of the recipe ablation."""
    main = {split: seed_means(results / 'main', split) for split in ('test', 'development')}
    main['test'].update({k: v for k, v in seed_means(results / 'explabels_training').items() if k[1] == 'hq_explabels'})
    m = lambda method, ts, metric, split='test': main[split].get((method, ts, metric))  # noqa: E731
    contrast = {**contrasts(results / 'explabels_training'), **contrasts(results / 'main')}

    # Supplementary Table 2: composition of the training sets and the benchmark
    sizes = {r['dataset']: r for r in read(results / 'dataset_sizes.tsv')}
    lines = []
    for key, name in (('original_transfew', 'TransFew original data'), ('original_deepgose', 'DeepGO-SE original data'),
                      ('original_dpfunc', 'DPFunc original data'), ('goa', 'HQ, GOA proteins only'), ('hq', 'HQ'),
                      ('hq_random', 'HQ + random'), ('archipelago', 'ArchipelaGO'), ('development', 'Development set'),
                      ('test', 'Test set')):
        r = sizes[key]
        if key == 'development':
            lines.append(r'\midrule')
        p, n = int(r['proteins']), int(r['propagated_pairs'])
        lines.append(f"{name} & {p:,} & {n:,} & {n / p:.2f} & {int(r['terms_used']):,} \\\\")
    write(out, 'table_composition.tex', lines)

    # Supplementary Table 3: every model on every training set
    lines = []
    for i, method in enumerate(['finch', 'transfew', 'deepgose', 'dpfunc', 'embedding_knn', 'diamond_knn', 'naive']):
        if i:
            lines.append(r'\addlinespace')
        first = True
        for ts, name in SETS:
            a = m(method, ts, 'auprc')
            if a is None:
                continue
            cells = [f"{r3(a['value'])} {ci(a['ci_low'], a['ci_high'])}"]
            cells += [r3(m(method, ts, k)['value']) for k in ('fmax', 'f1', 'protein_fmax')]
            cells.append(f"{float(m(method, ts, 'smin')['value']):.2f}")
            lines.append(f"{NAMES[method] if first else ''} & {name} & " + ' & '.join(cells) + r' \\')
            first = False
    write(out, 'table_ladder.tex', lines)

    # Supplementary Table 4: paired differences
    def differences(table, method, a, method_b, b):
        cells = []
        for k, places in (('auprc', 3), ('fmax', 3), ('smin', 2)):
            r = table[k, method, a, method_b, b]
            cells.append(f"{r3(r['difference'], True, places)} {ci(r['ci_low'], r['ci_high'], places)}")
        return cells

    lines = []
    for method in PUBLISHED:
        for a, b in (('hq_explabels', 'original'), ('hq', 'hq_explabels'), ('hq', 'goa'), ('hq', 'original'),
                     ('archipelago', 'hq'), ('archipelago', 'original')):
            lines.append(f'{NAMES[method]} & {SETNAME[a]} $-$ {SETNAME[b]} & '
                         + ' & '.join(differences(contrast, method, a, method, b)) + r' \\')
        lines.append(r'\addlinespace')
    for a, b in (('hq_random', 'hq'), ('archipelago', 'hq'), ('archipelago', 'hq_random')):
        lines.append(f'Finch & {SETNAME[a]} $-$ {SETNAME[b]} & '
                     + ' & '.join(differences(contrast, 'finch', a, 'finch', b)) + r' \\')
    lines.append(r'\addlinespace')
    for method in PUBLISHED:
        lines.append(f'Finch $-$ {NAMES[method]} & on ArchipelaGO & '
                     + ' & '.join(differences(contrast, 'finch', 'archipelago', method, 'archipelago')) + r' \\')
    write(out, 'table_contrasts.tex', lines)

    # TransFew's checkpoint rule: last epoch (as published), lowest development BCE, highest development AUPRC
    lines = []
    variants = ('transfew', 'transfew_bcebest', 'transfew_auprcbest')
    for ts, name in SETS:
        values = [m(v, ts, 'auprc') for v in variants]
        if None in values:
            continue
        devs = [m(v, ts, 'auprc', 'development') for v in variants]
        lines.append(f'{name} & ' + ' & '.join(f"{r3(v['value'])} {ci(v['ci_low'], v['ci_high'])}" for v in values)
                     + ' & ' + ' & '.join(r3(v['value']) for v in devs) + r' \\')
    write(out, 'table_transfew.tex', lines)

    label_counts(results, out)

    # test proteins without training relatives
    folder = results / 'no_hit'
    subsets = [('none', 5000, 2734)] + [(f"{float(r['cutoff']):g}", int(r['proteins']), int(r['clusters']))
                                         for r in read(folder / 'subsets.tsv')]
    nohit = {(r['subset'], r['method'], r['base']): r for r in read(folder / 'differences.tsv')}
    cutname = {'none': 'all test proteins', '1e-10': r'$10^{-10}$', '1e-06': r'$10^{-6}$', '0.001': r'$10^{-3}$',
               '0.01': r'$10^{-2}$', '0.1': r'$10^{-1}$', '1': '1', '10': '10'}
    lines = []
    for cut, n, k in subsets:
        random_gain = (contrast if cut == 'none' else contrasts(folder / f'E{cut}'))['auprc', 'finch', 'hq_random', 'finch', 'hq']
        rows = [nohit[cut, 'finch', 'hq'], nohit[cut, 'finch', 'hq_random'], random_gain,
                *(nohit[cut, p, 'original'] for p in PUBLISHED)]
        cells = [f"{r3(r['difference'], True)} {ci(r['ci_low'], r['ci_high'])}" for r in rows]
        lines.append(f'{cutname[cut]} & {n:,} & {k:,} & ' + ' & '.join(cells) + r' \\')
    write(out, 'table_nohit.tex', lines)

    # selection rules and number of picks (Finch)
    selection_contrast = contrasts(results / 'selection')
    lines = []
    for (a, b), name in ((('hq_random', 'hq'), 'HQ + random $-$ HQ'), (('pool_random', 'hq'), 'Random representatives $-$ HQ'),
                         (('filter_random', 'hq'), 'Filtered random representatives $-$ HQ'),
                         (('coverage_unfiltered', 'hq'), 'Ranking alone $-$ HQ'), (('archipelago', 'hq'), 'ArchipelaGO $-$ HQ'),
                         (('pool_random', 'hq_random'), 'One per cluster (random representatives $-$ HQ + random)'),
                         (('filter_random', 'pool_random'),
                          'Rare-function filter (filtered $-$ unfiltered random representatives)'),
                         (('archipelago', 'filter_random'),
                          'Coverage ranking (ArchipelaGO $-$ filtered random representatives)'),
                         (('archipelago_300k', 'archipelago'), '300,000 $-$ 150,000 picks'),
                         (('archipelago_600k', 'archipelago'), '600,000 $-$ 150,000 picks')):
        if (a, b) in (('pool_random', 'hq_random'), ('archipelago_300k', 'archipelago')):
            lines.append(r'\midrule')
        lines.append(f'{name} & ' + ' & '.join(differences(selection_contrast, 'finch', a, 'finch', b)) + r' \\')
    write(out, 'table_arms_diff.tex', lines)


def label_counts(results, out):
    """The counts of the experimental-only test labels (explabels_counts.tex)."""
    counts = {r['variant']: r for r in read(results / 'experimental_test_labels' / 'counts.tsv')}
    with open(out / 'explabels_counts.tex', 'w', encoding='utf-8', newline='\n') as handle:
        for variant, key in (('all', 'All'), ('exp', 'Exp')):
            r = counts[variant]
            handle.write(f"\\newcommand{{\\Count{key}Proteins}}{{{int(r['proteins']):,}}}\n"
                         f"\\newcommand{{\\Count{key}Pairs}}{{{int(r['propagated_pairs']):,}}}\n"
                         f"\\newcommand{{\\Count{key}Clusters}}{{{int(r['clusters']):,}}}\n"
                         f"\\newcommand{{\\Count{key}Terms}}{{{int(r['terms_with_positive']):,}}}\n")
        handle.write(f"\\newcommand{{\\CountExpDropped}}{{{int(counts['exp']['proteins_dropped']):,}}}\n")


def cost_table(results, out):
    """Training cost (table_cost.tex, from retrain.scoring.training_cost): per method, GPU-hours of one run on ArchipelaGO
    (the mean over its runs), of all its runs on ArchipelaGO, and of all its runs on every training set of the main
    comparison; methods from the cheapest run."""
    rows = read(results / 'training_cost.tsv')
    archipelago = {r['method']: r for r in rows if r['training_set'] == 'archipelago'}
    one = {m: float(r['gpu_hours']) / int(r['runs']) for m, r in archipelago.items()}
    lines = []
    for method in sorted(one, key=one.get):
        total = sum(float(r['gpu_hours']) for r in rows if r['method'] == method)
        lines.append(f"{NAMES[method]} & {half_up(one[method], places=2 if one[method] < 1 else 1)} & "
                     f"{half_up(archipelago[method]['gpu_hours'], places=1)} & {half_up(total, places=1)} \\\\")
    write(out, 'table_cost.tex', lines)


def recipe_tables(results, out):
    """The Finch recipe ablation (Fig. 4a) and the other models the main text compares, all trained on ArchipelaGO: test
    AUPRC, Fmax and Smin, development AUPRC, for the ablation the paired gain over the step it builds on, and trainable
    parameters in millions; then the paired comparisons cited in the main text."""
    folder = results / 'ladder'
    ladder = {split: seed_means(folder, split) for split in ('test', 'development')}
    ladder_contrast = contrasts(folder)
    settings = finch_config()
    models = {**settings['ladder'], **settings['variants']}
    parameters = json.loads((results / 'parameters.json').read_text(encoding='utf-8'))
    millions = lambda n: f'{n / 1e6:,.1f}'  # noqa: E731

    def setting(readout):
        """A readout's learning rate and dropout, in LaTeX."""
        r = settings['readouts'][readout]
        mantissa, exponent = f"{r['lr']:.0e}".split('e')
        lr = r'$10^{%d}$' % int(exponent) if mantissa == '1' else r'$%s\times10^{%d}$' % (mantissa, int(exponent))
        return lr + ('' if r['dropout'] is None else f", {r['dropout']:g}")

    def scores(model):
        t = lambda k: ladder['test'][model, 'archipelago', k]  # noqa: E731
        return (f"{r3(t('auprc')['value'])} {ci(t('auprc')['ci_low'], t('auprc')['ci_high'])} & "
                f"{r3(t('fmax')['value'])} & {float(t('smin')['value']):.2f} & "
                f"{r3(ladder['development'][model, 'archipelago', 'auprc']['value'])}")

    # where None: the learning rate and dropout of the model's one readout
    lines = []
    for model, name, base, where in (
            ('ladder_linear_esm2_3b', 'Linear readout on ESM2-3B', None, None),
            ('ladder_linear_esm2_15b', 'Linear readout on ESM2-15B', 'ladder_linear_esm2_3b', None),
            ('ladder_narrow', '+ hidden layer, 1,024 units, light dropout', 'ladder_linear_esm2_15b', None),
            ('ladder_hidden', '+ width 8,192 and heavy dropout', 'ladder_narrow', None),
            ('ladder_asam', '+ ASAM', 'ladder_hidden', None),
            ('finch', '+ ProtT5 network, averaged (Finch)', 'ladder_asam', 'Methods'),
            ('ladder_two_esm2_15b', 'Control: + a second ESM2-15B network, averaged', 'ladder_asam', 'as + ASAM')):
        if model == 'ladder_two_esm2_15b':
            lines.append(r'\midrule')
        step = DASH
        if base:
            d = ladder_contrast['auprc', model, 'archipelago', base, 'archipelago']
            step = f"{r3(d['difference'], True)} {ci(d['ci_low'], d['ci_high'])}"
        lines.append(f"{name} & {where or setting(next(iter(models[model])))} & {scores(model)} & {step} & "
                     f"{millions(sum(parameters[model].values()))} \\\\")
    lines.append(r'\midrule')
    for model, name in (('ladder_hidden_esm2_3b', 'Wide hidden layer on ESM2-3B'), ('ladder_hidden_esm1b', 'Wide hidden layer on ESM-1b')):
        lines.append(f"{name} & {setting(next(iter(models[model])))} & {scores(model)} & {DASH} & "
                     f"{millions(sum(parameters[model].values()))} \\\\")
    lines.append(r'\midrule')
    for model in ('deepgose', 'transfew', 'dpfunc'):
        parts = parameters[model].values()
        count = f'{millions(sum(parts) / len(parts))} per model' if model == 'deepgose' else millions(sum(parts))  # one part per world
        lines.append(f"{NAMES[model]} & published & {scores(model)} & {DASH} & {count} \\\\")
    write(out, 'table_recipe.tex', lines)

    # the paired comparisons cited in the main text
    lines = []
    for a, b, name in (('ladder_linear_esm2_15b', 'deepgose', 'Linear readout on ESM2-15B $-$ DeepGO-SE'),
                       ('ladder_narrow', 'deepgose', 'Hidden layer of 1,024 units on ESM2-15B $-$ DeepGO-SE'),
                       ('ladder_narrow', 'transfew', 'Hidden layer of 1,024 units on ESM2-15B $-$ TransFew'),
                       ('ladder_hidden_esm2_3b', 'deepgose', 'Wide hidden layer on ESM2-3B $-$ DeepGO-SE'),
                       ('ladder_hidden', 'transfew', 'Wide hidden layer on ESM2-15B $-$ TransFew'),
                       ('ladder_hidden_esm1b', 'dpfunc', 'Wide hidden layer on ESM-1b $-$ DPFunc'),
                       ('ladder_linear_esm2_15b', 'dpfunc', 'Linear readout on ESM2-15B $-$ DPFunc'),
                       ('finch', 'finch_prott5_l24', 'Finch $-$ Finch with ProtT5 layer 24')):
        cells = [paired(ladder_contrast[k, a, 'archipelago', b, 'archipelago'], places)
                 for k, places in (('auprc', 3), ('fmax', 3), ('smin', 2))]
        lines.append(f'{name} & ' + ' & '.join(cells) + r' \\')
    write(out, 'table_pairs.tex', lines)


def shared(results, out):
    """Finch and the published methods on every training set they share: table_shared.tex, and the precision-recall curves
    by training set (si_fig_pr_sets) and by model (si_fig_pr_models), each dot a curve's Fmax point."""
    main, trained = seed_means(results / 'main'), seed_means(results / 'explabels_training')
    main_pairs, trained_pairs = contrasts(results / 'main'), contrasts(results / 'explabels_training')
    value = lambda metric, m, s, k: float(metric[m, s, k]['value'])  # noqa: E731
    lines = []
    for s, title, own in SHARED:
        metric, pairs = (trained, trained_pairs) if s == 'hq_explabels' else (main, main_pairs)
        ref, ref_set = (own, 'original') if own else (max(PUBLISHED, key=lambda p: value(metric, p, s, 'auprc')), s)
        for k, label, places in (('auprc', 'Micro-AUPRC', 3), ('fmax', 'Fmax', 3), ('smin', 'Smin', 2)):
            lines.append(f"{title if k == 'auprc' else ''} & {NAMES[ref] if k == 'auprc' else ''} & {label} & "
                         f"{r3(value(metric, 'finch', s, k), places=places)} & {r3(value(metric, ref, ref_set, k), places=places)} & "
                         f"{paired(pairs[k, 'finch', s, ref, ref_set], places)}" + r' \\')
        lines.append(r'\addlinespace')
    write(out, 'table_shared.tex', lines[:-1])

    curves = {key: np.array(v).T for key, v in pr_curves(results).items()}

    def fmax_point(ax, recall, precision, color):
        k = int(np.argmax(2 * recall * precision / np.maximum(recall + precision, 1e-12)))
        ax.scatter(recall[k], precision[k], s=9, color=color, edgecolors=INK, linewidths=0.4, zorder=5)

    def square(ax):
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.02)
        ax.set_xticks(np.arange(0, 1.01, 0.2))
        ax.set_yticks(np.arange(0, 1.01, 0.2))
        ax.tick_params(labelsize=fs(6))
        ax.set_box_aspect(1)

    fig, grid = plt.subplots(2, 4, figsize=(DOUBLE, 3.75), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.065, right=0.995, bottom=0.095, top=0.955, wspace=0.08, hspace=0.2)
    grid[0, 3].set_visible(False)  # the top row holds the three original sets, the bottom row the sets all methods share
    for ax, (s, title, own) in zip([grid[0, 0], grid[0, 1], grid[0, 2], *grid[1]], SHARED):
        metric = trained if s == 'hq_explabels' else main
        series = [('finch', s, COLOURS['archipelago'], 1.4)] + \
                 [(p, 'original' if own else s, SHADES[p], 1.0) for p in ([own] if own else PUBLISHED)]
        for m, training_set, color, lw in series:
            recall, precision = curves[m, training_set]
            ax.plot(recall, precision, color=color, lw=lw, solid_capstyle='butt', zorder=4 if m == 'finch' else 3,
                    label=f"{NAMES[m]} {r3(value(metric, m, training_set, 'auprc'))} | {r3(value(metric, m, training_set, 'fmax'))}")
            fmax_point(ax, recall, precision, color)
        ax.legend(loc='lower left', title='AUPRC | Fmax', frameon=False, fontsize=fs(5.5), title_fontsize=fs(5.5),
                  handlelength=0.9, handletextpad=0.3, alignment='left', borderaxespad=0.1, labelspacing=0.2)
        ax.set_title(title, fontsize=fs(6.5), fontweight='normal', loc='left', pad=2.5)
        square(ax)
        ax.tick_params(axis='x', labelbottom=True)
    for ax in (grid[0, 0], grid[1, 0]):
        ax.set_ylabel('Precision', fontsize=fs(7))
    for ax in grid[1]:
        ax.set_xlabel('Recall', fontsize=fs(7))
    save(fig, out, 'si_fig_pr_sets')

    fig, grid = plt.subplots(1, 4, figsize=(DOUBLE, 2.15), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.065, right=0.995, bottom=0.17, top=0.8, wspace=0.08)
    for ax, m in zip(grid, ('finch', 'deepgose', 'transfew', 'dpfunc')):
        for s, title, own in SHARED:
            if m == 'finch' or own in (None, m):
                recall, precision = curves[m, 'original' if own and m != 'finch' else s]
                ax.plot(recall, precision, color=SET_LINE[s], lw=1.4 if s == 'archipelago' else 1.0, solid_capstyle='butt',
                        label=title)
                fmax_point(ax, recall, precision, SET_LINE[s])
        ax.set_title(NAMES[m], fontsize=fs(6.5), fontweight='normal', loc='left', pad=2.5)
        ax.set_xlabel('Recall', fontsize=fs(7))
        square(ax)
    grid[0].set_ylabel('Precision', fontsize=fs(7))
    handles, labels = grid[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=4, frameon=False, fontsize=fs(6), handlelength=1.6, columnspacing=1.2,
               labelspacing=0.25, borderaxespad=0.1)
    save(fig, out, 'si_fig_pr_models')


def draw(results, out):
    out = out / 'si'
    out.mkdir(parents=True, exist_ok=True)
    tables(results, out)
    recipe_tables(results, out)
    cost_table(results, out)
    shared(results, out)
