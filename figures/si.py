"""Tables and figures of the Supplementary Information: LaTeX table bodies (rows only; the SI typesets the headers) and
the figures, written into <out>/si/. Reads, under <root>/results/:

  main/{metrics.tsv,contrasts.tsv}                         every model on every training set (retrain.scoring.score main)
  dataset_sizes.tsv                                        retrain.scoring.sizes
  experimental_test_labels/counts.tsv                      retrain.scoring.test_labels experimental
  explabels_training/{metrics.tsv,contrasts.tsv}           retrain.scoring.score explabels_training
  no_hit/{subsets.tsv,differences.tsv,E<cutoff>/}          retrain.scoring.no_hit
  pr_curves/curves.tsv                                     retrain.scoring.pr_curves
  selection/{metrics.tsv,contrasts.tsv}                    retrain.scoring.score selection
  term_strata/gains.tsv                                    retrain.scoring.term_strata
  ladder/{metrics.tsv,contrasts.tsv}                       retrain.scoring.score ladder
  parameters.json                                          retrain.scoring.parameters
  training_cost.tsv                                        retrain.scoring.training_cost
  macro/contrasts.tsv                                      retrain.scoring.macro (Supplementary Table 15)
  outside_masked/<method>/{table,contrasts}.tsv            retrain.outside.masked bootstrap (its work folder; Supplementary
                                                           Table 17, Supplementary Figure 3)

Main-text Table 2 with intervals (table2_intervals.tex) is written beside Table 2 by figures/main.py.

Supplementary Table 1 (the annotation sources) is typed in the SI and the label-review prompt is released with the
reviewed labels; neither is generated here.
"""
import json

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.legend_handler import HandlerTuple
from matplotlib.transforms import blended_transform_factory

from figures.style import (COLOURS, DOUBLE, GREY, INK, NAMES, PUBLISHED, contrasts, fs, half_up, pr_curves, read, save,
                           seed_means)
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
# the published methods' own benchmarks (Supplementary Note 14): (method, test set, name)
BENCHMARKS = (('deepgose', 'test', 'DeepGO-SE'), ('transfew', 'Test_all', r'TransFew, Test\_all'), ('dpfunc', 'test', 'DPFunc'))
ORDER = ('deepgose', 'transfew', 'dpfunc')  # main-text Table 2's order
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


def sd(value):
    """A run-to-run standard deviation; models with one fit or one ensemble have none."""
    return DASH if value in ('', None) else f'{float(value):.4f}'


def write(out, name, lines):
    with open(out / name, 'w', encoding='utf-8', newline='\n') as handle:
        handle.write('\n'.join(lines) + '\n')


def tables(results, out):
    """Every table body of the SI but those of the recipe ablation."""
    main = {split: seed_means(results / 'main', split) for split in ('test', 'development')}
    # not DIAMOND-kNN: that run scores it at the screen's cutoff only, and the paper reports it without a cutoff
    main['test'].update({k: v for k, v in seed_means(results / 'explabels_training').items()
                         if k[1] == 'hq_explabels' and k[0] != 'diamond_knn'})
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
            cells = [f"{r3(a['value'])} {ci(a['ci_low'], a['ci_high'])}", sd(a['seed_sd'])]
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
    runs = {}
    for r in read(results / 'selection' / 'metrics.tsv'):
        if r['method'] == 'finch' and r['split'] == 'test' and r['metric'] == 'auprc' and r['seed'].isdigit():
            runs.setdefault(r['training_set'], {})[r['seed']] = float(r['value'])
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
        cells = differences(selection_contrast, 'finch', a, 'finch', b)
        cells.insert(1, sd(np.std([runs[a][s] - runs[b][s] for s in sorted(runs[a])], ddof=1)))  # the same seed in both
        lines.append(f'{name} & ' + ' & '.join(cells) + r' \\')
    write(out, 'table_arms_diff.tex', lines)

    # the main gains on common and rare terms
    strata = {(r['stratum'], r['method'], r['training_set'], r['minus_training_set']): r
              for r in read(results / 'term_strata' / 'gains.tsv')}
    lines = []
    for method, a, b in [(p, 'archipelago', 'original') for p in PUBLISHED] + [('finch', 'archipelago', 'hq_random')]:
        lines.append(f'{NAMES[method]} & {SETNAME[a]} $-$ {SETNAME[b]} & '
                     + ' & '.join(paired(strata[s, method, a, b], 3) for s in ('all', 'common', 'rare')) + r' \\')
    write(out, 'table_term_strata.tex', lines)


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
        return (f"{r3(t('auprc')['value'])} {ci(t('auprc')['ci_low'], t('auprc')['ci_high'])} & {sd(t('auprc')['seed_sd'])} & "
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


def macro_table(results, out):
    """Supplementary Table 15: the paper's main comparisons in micro-AUPRC (paired cluster bootstrap) and in term-centric
    AUPRC (paired bootstrap over terms), from retrain.scoring.macro (results/macro/contrasts.tsv; its cluster bootstrap of
    macro-AUPRC is a diagnostic and not read)."""
    rows = [r for r in read(results / 'macro' / 'contrasts.tsv') if r['bootstrap'] in ('clusters', 'terms')]

    def cells(pick):
        return ' & '.join(paired(next(r for r in rows if r['metric'] == k and pick(r)), 3) for k in ('micro_auprc', 'macro_auprc'))

    def pair(method, training_set, minus_method, minus_training_set):
        return lambda r: (r['method'], r['training_set'], r['minus_method'], r['minus_training_set']) == \
            (method, training_set, minus_method, minus_training_set)

    lines = [f"{NAMES[m]} & ArchipelaGO $-$ Original & {cells(pair(m, 'archipelago', m, 'original'))}" + r' \\' for m in ORDER]
    lines += [r'\addlinespace',
              f"Finch & ArchipelaGO $-$ HQ + random & {cells(pair('finch', 'archipelago', 'finch', 'hq_random'))}" + r' \\',
              r'\addlinespace']
    lines += [f"Finch $-$ {NAMES[m]} & on ArchipelaGO & {cells(pair('finch', 'archipelago', m, 'archipelago'))}" + r' \\'
              for m in ORDER]
    # the spread: the best minus the worst published method on ArchipelaGO, by each metric's own point estimates
    lines += [r'\addlinespace', 'Spread & Best $-$ worst method, on ArchipelaGO & '
              + cells(lambda r: r['comparison'].startswith('Spread')) + r' \\']
    write(out, 'table_macro.tex', lines)


# Supplementary Note 14: Finch on the published methods' benchmarks (retrain.outside.masked's work folder,
# results/outside_masked/<method>/{table,contrasts}.tsv): micro-AUPRC (seed means) on all test proteins of each benchmark, on
# the paper's own test labels (unmasked) and with the test proteins' computational GOA annotations masked, with paired 95%
# intervals from a bootstrap over the test sequences. Models: method, the method retrained on its training split; finch_A,
# Finch on that split; finch_A_screened, Finch on that split without the test proteins' relatives (TransFew and DPFunc);
# finch_B, Finch on ArchipelaGO screened against the test proteins. TransFew: its Test_all.
VERSIONS = (('unmasked', "Paper's"), ('computational', 'Masked'))


def outside(folder, method, test_set, version):
    """{model: table.tsv row}, {contrast: contrasts.tsv row of micro-AUPRC} of all test proteins and one label version."""
    cell = {r['model']: r for r in read(folder / method / 'table.tsv')
            if (r['test_set'], r['subset'], r['version']) == (test_set, 'all', version)}
    diff = {r['contrast']: r for r in read(folder / method / 'contrasts.tsv')
            if (r['test_set'], r['subset'], r['version'], r['metric']) == (test_set, 'all', version, 'micro_AUPRC')}
    return cell, diff


def with_interval(row):
    """A paired difference followed by its interval in small grey type (the SI's \\ci)."""
    value, interval = paired(row, 3).split(' ', 1)
    return rf'{value}\ci{{{interval}}}'


def block(name, proteins, body):
    """One benchmark's rows: its name in the first row, its number of test proteins in the second."""
    first = [name, rf'\textcolor{{muted}}{{{proteins:,} proteins}}'] + [''] * (len(body) - 2)
    return [' & '.join([b, *r]) + r' \\' for b, r in zip(first, body)]


# Supplementary Table 17b: Finch on ArchipelaGO against Finch on the method's training split, as released and without the
# test proteins' relatives (DeepGO-SE: no screened run)
SPLITS = (('finch_A', 'As released', 'finch_B - finch_A'),
          ('finch_A_screened', 'Relatives removed', 'finch_B - finch_A_screened'))


def outside_table(results, out):
    """Supplementary Table 17: Finch on each published method's own benchmark (table_outside_methods.tex, a: Finch and the
    method retrained, both trained on the method's training split; table_outside_data.tex, b: Finch on ArchipelaGO against
    Finch on that split)."""
    folder = results / 'outside_masked'
    lines = []
    for method, test_set, name in BENCHMARKS:
        body = []
        for version, label in VERSIONS:
            cell, diff = outside(folder, method, test_set, version)
            body.append([label, r3(cell['method']['micro_AUPRC_mean']), r3(cell['finch_A']['micro_AUPRC_mean']),
                         with_interval(diff['finch_A - method'])])
        lines += ([r'\addlinespace'] if lines else []) + block(name, int(cell['method']['proteins']), body)
    write(out, 'table_outside_methods.tex', lines)

    lines = []
    for method, test_set, name in BENCHMARKS:
        body = []
        for version, label in VERSIONS:
            cell, diff = outside(folder, method, test_set, version)
            present = [s for s in SPLITS if s[0] in cell]
            body += [[label if k == 0 else '', split, r3(cell[model]['micro_AUPRC_mean']), r3(cell['finch_B']['micro_AUPRC_mean']),
                      with_interval(diff[contrast])] for k, (model, split, contrast) in enumerate(present)]
        lines += ([r'\addlinespace'] if lines else []) + block(name, int(cell['method']['proteins']), body)
    write(out, 'table_outside_data.tex', lines)


# Supplementary Figure 3, in the look of main-text Fig. 4b (figures/main.py). a: Finch against the method retrained, both
# trained on the method's training split, on the paper's labels, Finch's lead at the right. b: Finch on ArchipelaGO minus
# Finch on the method's data, with the paired interval, on the paper's labels, with the computational annotations masked,
# and masked against the method's data without the test proteins' relatives (TransFew and DPFunc). One colour, Finch's
# violet; the fill counts the corrections. Drawn at the SI's text width, so that 7 pt is 7 pt on the page.
MARK = {'deepgose': ('o', 4.2), 'transfew': ('s', 3.8), 'dpfunc': ('^', 4.6)}  # as Fig. 4b
FINCH, DOT_INK, ROW_LINE, ZERO_LINE, GRID = COLOURS['archipelago'], '#55595f', '#e6e7e9', '#8f939a', '#cfd1d5'
TEXT = 7  # points, every text of the figure
WIDTH, LEFT, RIGHT = 6.85, 1.25, 0.55  # inches: the SI's text width (174 mm); the row names; panel a's leads
STEPS = (('unmasked', 'finch_B - finch_A', "Paper's labels", dict(mfc='white')),
         ('computational', 'finch_B - finch_A', 'Computational annotations masked', dict(fillstyle='left', mfc=FINCH,
                                                                                         markerfacecoloralt='white')),
         ('computational', 'finch_B - finch_A_screened', 'Masked, and relatives removed from own data', dict(mfc=FINCH)))
SUBROW = (0.24, 0.0, -0.24)  # panel b: each step's place within its row
RULE = '#b9bcc2'  # panel b: the rules between the benchmarks
ROW_A, ROW_B = 0.3, 0.44  # inches per row
BOTTOM, KEY, GAP = 0.42, 0.24, 0.1  # inches: x axis and its label; the key above each panel; between the panels


def signed_text(value):
    return half_up(value, sign=True).replace('-', '−')


def rows(fig, ax, folder):
    """Each benchmark's row: its name and number of test proteins at the figure's left, top row first."""
    left = blended_transform_factory(fig.transFigure, ax.transData)
    for i, (method, test_set, name) in enumerate(BENCHMARKS):
        cell, _ = outside(folder, method, test_set, 'unmasked')
        for text, dy, color in ((name.replace(r'\_', '_'), 4.2, 'black'), (f"{int(cell['method']['proteins']):,} test proteins", -4.2, GREY)):
            ax.annotate(text, (0.06 / WIDTH, 2 - i), xycoords=left, xytext=(0, dy), textcoords='offset points', ha='left',
                        va='center', fontsize=TEXT, color=color)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_ylim(-0.6, 2.6)
    ax.tick_params(axis='x', labelsize=TEXT)


def ticks(ax, values, label):
    """The x ticks, each with a light dashed line up through the panel."""
    values = np.round(values, 2)
    ax.set_xticks(values, [label(v) for v in values])
    ax.grid(axis='x', color=GRID, ls=(0, (3, 2)), lw=0.5)
    ax.set_axisbelow(True)


def key(fig, ax, letter, entries):
    """The panel letter at the figure's left and the panel's key above its axes, in one row."""
    handle = lambda s: plt.Line2D([], [], **s)  # noqa: E731
    ax.legend([tuple(handle(s) for s in style) if isinstance(style, list) else handle(style) for _, style in entries],
              [label for label, _ in entries], loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=len(entries), frameon=False,
              fontsize=TEXT, handletextpad=0.4, columnspacing=1.6, borderaxespad=0.25, borderpad=0,
              handler_map={tuple: HandlerTuple(ndivide=None, pad=0.2)})
    fig.text(0.06 / WIDTH, ax.get_position().y1 + 0.1 / fig.get_figheight(), letter, fontsize=8, fontweight='bold',
             ha='left', va='center')


def outside_figure(results, out):
    """Supplementary Figure 3 (si_fig_outside)."""
    folder = results / 'outside_masked'
    height_b, height_a = 3.2 * ROW_B, 3.2 * ROW_A
    height = BOTTOM + height_b + KEY + GAP + BOTTOM + height_a + KEY
    fig = plt.figure(figsize=(WIDTH, height))
    inch = lambda x, y, w, h: [x / WIDTH, y / height, w / WIDTH, h / height]  # noqa: E731
    a = fig.add_axes(inch(LEFT, height - KEY - height_a, WIDTH - LEFT - RIGHT, height_a))
    b = fig.add_axes(inch(LEFT, BOTTOM, WIDTH - LEFT - RIGHT, height_b))

    # a: Finch and the method, both on the method's data
    rows(fig, a, folder)
    for i, (method, test_set, _) in enumerate(BENCHMARKS):
        y = 2 - i
        cell, diff = outside(folder, method, test_set, 'unmasked')
        own, finch = float(cell['method']['micro_AUPRC_mean']), float(cell['finch_A']['micro_AUPRC_mean'])
        a.axhline(y, color=ROW_LINE, lw=0.6, zorder=1)
        a.plot([own, finch], [y, y], color=FINCH, lw=1.4, alpha=0.35, solid_capstyle='butt', zorder=2)
        a.plot(own, y, MARK[method][0], ms=MARK[method][1], mfc=DOT_INK, mec='white', mew=0.5, zorder=3)
        a.plot(finch, y, 'o', ms=5.6, mfc=FINCH, mec='white', mew=0.5, zorder=4)
        a.text(finch, y + 0.24, r3(finch), ha='center', va='bottom', fontsize=TEXT, zorder=5)
        a.text(1.015, y, signed_text(diff['finch_A - method']['difference']), transform=a.get_yaxis_transform(), ha='left',
               va='center', fontsize=TEXT)
    a.annotate('Lead', (1.015, 1.0), xycoords='axes fraction', xytext=(0, 5.5), textcoords='offset points', ha='left',
               va='center', fontsize=TEXT, color=GREY)  # level with the key
    a.set_xlim(0.43, 0.67)
    ticks(a, np.arange(0.44, 0.661, 0.02), lambda v: f'{v:.2f}')
    a.set_xlabel("Micro-AUPRC, paper's test labels", fontsize=TEXT)
    key(fig, a, 'a', [('Finch', dict(ls='', marker='o', ms=5.6, mfc=FINCH, mec='white', mew=0.5)),
                      ("The method, retrained by us (both trained on the method's own data)",
                       [dict(ls='', marker=MARK[m][0], ms=MARK[m][1], mfc=DOT_INK, mec='white', mew=0.5) for m in ORDER])])

    # b: what closes ArchipelaGO's gap
    rows(fig, b, folder)
    for y in (1.5, 0.5):  # a rule between the benchmarks, from the row names to the panel's right edge
        fig.add_artist(plt.Line2D([0.06 / WIDTH, (WIDTH - RIGHT) / WIDTH], [y, y], color=RULE, lw=0.6,
                                  transform=blended_transform_factory(fig.transFigure, b.transData)))
    b.axvline(0, color=ZERO_LINE, lw=0.7, zorder=1)
    for i, (method, test_set, _) in enumerate(BENCHMARKS):
        for (version, contrast, _, style), dy in zip(STEPS, SUBROW):
            _, diff = outside(folder, method, test_set, version)
            if contrast not in diff:  # DeepGO-SE: no screened run
                continue
            d, y = diff[contrast], 2 - i + dy
            b.plot([float(d['ci_low']), float(d['ci_high'])], [y, y], color=FINCH, lw=1.0, solid_capstyle='butt', zorder=2)
            b.plot(float(d['difference']), y, 'o', ms=5.6, mec=FINCH, mew=1.0, zorder=3, **style)
    b.set_xlim(-0.27, 0.17)
    ticks(b, np.arange(-0.25, 0.151, 0.05), lambda v: '0' if v == 0 else f'{v:+.2f}'.replace('-', '−'))
    b.set_xlabel("Finch trained on ArchipelaGO minus Finch trained on the method's own data (micro-AUPRC)", fontsize=TEXT)
    key(fig, b, 'b', [(label, dict(color=FINCH, lw=1.0, marker='o', ms=5.6, mec=FINCH, mew=1.0, **style))
                      for _, _, label, style in STEPS])
    save(fig, out, 'si_fig_outside')


def draw(results, out):
    out = out / 'si'
    out.mkdir(parents=True, exist_ok=True)
    tables(results, out)
    recipe_tables(results, out)
    cost_table(results, out)
    shared(results, out)
    macro_table(results, out)
    outside_table(results, out)
    outside_figure(results, out)
