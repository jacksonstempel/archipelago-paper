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
  outside_masked/<method>/{table,contrasts}.tsv            retrain.outside.masked bootstrap (its work folder; Supplementary
                                                           Tables 16-17, Figures 3-4)
  outside_masked/deepgose/pooling_contrasts.tsv            retrain.outside.masked pooling deepgose (Supplementary Table 19)
  outside_compare/<method>/relatives.tsv                   retrain.outside.compare (its outside/compare/relatives.tsv;
                                                           Supplementary Table 18)
  shared_terms/summary.tsv                                 retrain.scoring.shared_terms (Supplementary Table 20)

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
# the published methods' own benchmarks (Supplementary Note 14): (method, name, test set)
BENCHMARKS = (('deepgose', 'DeepGO-SE', 'test'), ('transfew', r'TransFew, Test\_all', 'Test_all'),
              ('transfew', r'TransFew, Test\_novel', 'Test_novel'), ('dpfunc', 'DPFunc', 'test'))
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


def outside_cell(folder, method, test_set, subset, version, metric):
    """{model: table.tsv row}, {contrast: contrasts.tsv row} of one row of retrain.outside.masked's tables."""
    table, pairs = read(folder / method / 'table.tsv'), read(folder / method / 'contrasts.tsv')
    cell = {r['model']: r for r in table if (r['test_set'], r['subset'], r['version']) == (test_set, subset, version)}
    diff = {r['contrast']: r for r in pairs if (r['test_set'], r['subset'], r['version'], r['metric']) == (test_set, subset, version, metric)}
    return cell, diff


def outside_tables(results, out):
    """Supplementary Tables 16-19, Finch on the published methods' own benchmarks: micro-AUPRC and micro-Fmax by test
    labels and test proteins, with the paired differences (retrain.outside.masked's work folder, results/outside_masked);
    micro-AUPRC by each test protein's relation to the method's training split (retrain.outside.compare's relatives.tsv,
    results/outside_compare/<method>/); on DeepGO-SE's benchmark, Finch on ArchipelaGO minus Finch on DeepGO-SE's data by
    how each metric pools the pairs (retrain.outside.masked pooling deepgose)."""
    folder = results / 'outside_masked'
    models = ('finch_A', 'finch_A_screened', 'finch_B', 'method', 'released')  # A-screened: TransFew and DPFunc only
    for metric, name in (('micro_AUPRC', 'table_outside_auprc.tex'), ('micro_Fmax', 'table_outside_fmax.tex')):
        lines = []
        for method, label, test_set in BENCHMARKS:
            if lines:
                lines.append(r'\midrule')
            bench_lines = label.replace(', ', ',\n').split('\n')  # names over two rows, so that the table fits the page
            k = 0
            for version, vname in (('unmasked', "Paper's"), ('computational', 'Computational masked'), ('all_goa', 'All GOA masked')):
                version_lines = vname.replace(' masked', '\nmasked').split('\n')
                for j, (subset, sname) in enumerate((('all', 'All'), ('no_relative', 'No relative'))):
                    cell, diff = outside_cell(folder, method, test_set, subset, version, metric)
                    values = [r3(cell[m][metric + '_mean']) if m in cell else DASH for m in models]
                    screened = diff.get('finch_B - finch_A_screened')
                    lines.append(' & '.join([bench_lines[k] if k < len(bench_lines) else '',
                                             version_lines[j] if j < len(version_lines) else '',
                                             f"{sname} ({int(cell['finch_A']['proteins']):,})", *values,
                                             paired(diff['finch_B - finch_A'], 3), paired(screened, 3) if screened else DASH,
                                             paired(diff['finch_A - method'], 3)]) + r' \\')
                    k += 1
        write(out, name, lines)

    lines = []
    for method, label, test_set in BENCHMARKS:
        relatives = read(results / 'outside_compare' / method / 'relatives.tsv')
        if lines:
            lines.append(r'\midrule')
        first = True
        for group, name in (('identical', 'Identical sequence'), ('diamond', 'DIAMOND hit only'), ('none', 'No relative')):
            cell = {r['model']: r for r in relatives if (r['test_set'], r['group']) == (test_set, group)}
            if not cell:
                continue
            values = [r3(cell[m]['micro_AUPRC_mean']) if m in cell else DASH for m in ('finch_A', 'finch_B', 'method', 'released')]
            lines.append(' & '.join([label if first else '', f"{name} ({int(cell['finch_A']['proteins']):,})", *values]) + r' \\')
            first = False
    write(out, 'table_outside_relatives.tex', lines)

    pooling = {(r['version'], r['metric']): r for r in read(folder / 'deepgose' / 'pooling_contrasts.tsv')
               if r['test_set'] == 'test' and r['contrast'] == 'finch_B - finch_A'}
    lines = []
    for k, name in (('micro_AUPRC', 'Micro-AUPRC'), ('micro_Fmax', 'Micro-Fmax'), ('protein_Fmax', 'Protein-centric Fmax'),
                    ('macro_AUPRC', 'Term-centric AUPRC')):
        lines.append(f'{name} & ' + ' & '.join(paired(pooling[version, k], 3) for version in ('unmasked', 'computational', 'all_goa'))
                     + r' \\')
    write(out, 'table_outside_deepgose.tex', lines)


def vocabulary_table(results, out):
    """Supplementary Table 20: the data comparisons over every term and over the terms each original training set covers
    (retrain.scoring.shared_terms, results/shared_terms/summary.tsv)."""
    vocab = {(r['comparison'], r['terms'], r['metric']): r for r in read(results / 'shared_terms' / 'summary.tsv')}
    lines = []
    for method in ('deepgose', 'transfew', 'dpfunc'):
        own, name = f'original_{method}', NAMES[method]
        if lines:
            lines.append(r'\addlinespace')
        lines.append(rf"\multicolumn{{8}}{{@{{}}l}}{{\textit{{Terms of {name}'s original data: "
                     rf"{int(vocab[f'{name}, original to HQ', own, 'auprc']['n_terms']):,} of 8,739}}}} \\")
        for comparison, model, label in ((f'{name}, original to HQ', name, 'HQ $-$ Original'),
                                         (f'{name}, HQ to ArchipelaGO', '', 'ArchipelaGO $-$ HQ'),
                                         (f'{name}, original to ArchipelaGO', '', 'ArchipelaGO $-$ Original'),
                                         (f"Finch, {name}'s original set to ArchipelaGO", 'Finch', 'ArchipelaGO $-$ Original'),
                                         ('Spread of the published methods on ArchipelaGO', 'Spread', 'Best $-$ worst method')):
            cells = [paired(vocab[comparison, t, 'auprc'], 3) for t in ('all', own)]
            cells += [r3(vocab[comparison, t, 'fmax']['difference'], True) for t in ('all', own)]
            cells += [r3(vocab[comparison, t, 'smin']['difference'], True, 2) for t in ('all', own)]
            lines.append(f'{model} & {label} & ' + ' & '.join(cells) + r' \\')
    write(out, 'table_vocab.tex', lines)


# Supplementary Figures 3 and 4, in the look of main-text Fig. 4b: one row per benchmark and test labels, a dot per model,
# the differences printed to the right of each panel. Micro-AUPRC, seed means (retrain.outside.masked's tables).
MARK = {'deepgose': ('o', 4.2), 'transfew': ('s', 3.8), 'dpfunc': ('^', 4.6)}  # as Fig. 4b
DOT_INK, ROW_LINE, FINCH_MS = '#55595f', '#e6e7e9', 5.6
ROW_VERSIONS = (('unmasked', "Paper's labels"), ('computational', 'Computational masked'), ('all_goa', 'All GOA masked'))
GROUP_GAP, HEADER = 0.9, 0.8  # rows between benchmarks; the benchmark's name, above its first row
OUTSIDE_X = 0.25, 0.80
FINCH_OWN = dict(marker='o', ms=FINCH_MS, mfc=COLOURS['original'], mec='white', mew=0.5)
FINCH_SCREENED = dict(marker='o', ms=FINCH_MS - 0.8, mfc='white', mec=COLOURS['original'], mew=1.0)
FINCH_ARCHIPELAGO = dict(marker='o', ms=FINCH_MS, mfc=COLOURS['archipelago'], mec='white', mew=0.5)
FINCH = FINCH_ARCHIPELAGO  # Fig. 4b's Finch


def row_positions():
    """[(benchmark, version, version name, y)], top row first, and [(benchmark, y of its name)]: one row per test-label
    version, the benchmark's name above its rows."""
    rows, groups, y = [], [], 0.0
    for benchmark in BENCHMARKS:
        groups.append((benchmark, y + HEADER))
        for version, vname in ROW_VERSIONS:
            rows.append((benchmark, version, vname, y))
            y -= 1
        y -= GROUP_GAP
    return rows, groups


def dodge(values, width):
    """Vertical offsets for dots closer than `width` along x (Fig. 4b's rule for two near dots, for any number): each run
    of near dots spreads over the row, the higher value higher."""
    order = sorted(values, key=values.get)
    runs = [[order[0]]]
    for k in order[1:]:
        if values[k] - values[runs[-1][-1]] < width:
            runs[-1].append(k)
        else:
            runs.append([k])
    out = {}
    for run in runs:
        step = 0.34 if len(run) == 2 else 0.24  # two dots: +-0.17 of a row, as Fig. 4b
        for i, k in enumerate(run):
            out[k] = (i - (len(run) - 1) / 2) * step
    return out


def dot_panel(fig, ax, folder, rows, groups, subset, xlim, names=True):
    """Row lines, test-label names (names=True: and the benchmarks' names at the figure's left), each benchmark's number
    of test proteins in the subset, and the micro-AUPRC axis."""
    for *_, y in rows:
        ax.axhline(y, color=ROW_LINE, lw=0.6, zorder=1)
    ax.set_yticks([y for *_, y in rows], [v for _, _, v, _ in rows] if names else [], fontsize=fs(6.5))
    ax.tick_params(axis='y', length=0, pad=3)
    ax.spines['left'].set_visible(False)
    ax.set_ylim(rows[-1][3] - 0.6, rows[0][3] + HEADER + 0.45)
    ax.set_xlim(*xlim)
    ax.set_xticks(np.arange(np.ceil(xlim[0] * 10) / 10, xlim[1] + 0.001, 0.1))
    ax.tick_params(axis='x', labelsize=fs(6.5))
    ax.set_xlabel('Micro-AUPRC', fontsize=fs(7))
    left = blended_transform_factory(fig.transFigure, ax.transData)
    for (method, label, test_set), y in groups:
        if names:
            ax.text(0.06 / fig.get_figwidth(), y, label.replace(r'\_', '_'), transform=left, ha='left', va='center',
                    fontsize=fs(6.5))
        cell, _ = outside_cell(folder, method, test_set, subset, 'unmasked', 'micro_AUPRC')
        ax.text(0.005, y, f"{int(cell['finch_A']['proteins']):,} test proteins", transform=ax.get_yaxis_transform(),
                ha='left', va='center', fontsize=fs(6), color=GREY)


def dot_width(fig, ax, points=5.0):
    """`points` of the page as a length along the x axis."""
    lo, hi = ax.get_xlim()
    return points / 72 / (ax.get_position().width * fig.get_figwidth()) * (hi - lo)


def signed(row):
    return half_up(row['difference'], sign=True).replace('-', '\u2212')


def method_style(method, hollow=False):
    shape, ms = MARK[method]
    return dict(marker=shape, ms=ms - (0.5 if hollow else 0), mfc='white' if hollow else DOT_INK, mec=DOT_INK if hollow else 'white',
                mew=0.8 if hollow else 0.5)


def key(fig, entries):
    """The figure's key above the panels, in one row: (label, style or [styles]); a list shares one entry (one per method)."""
    handle = lambda s: plt.Line2D([], [], ls='', **s)  # noqa: E731
    fig.legend([tuple(handle(s) for s in style) if isinstance(style, list) else handle(style) for _, style in entries],
               [label for label, _ in entries], loc='upper left', bbox_to_anchor=(0.06 / fig.get_figwidth(), 0.995),
               ncol=len(entries), frameon=False, fontsize=fs(6.5), handletextpad=0.3, columnspacing=1.4, borderaxespad=0,
               borderpad=0, handler_map={tuple: HandlerTuple(ndivide=None, pad=0.15)})


def column_head(ax, x, text, **kw):
    """A head over a number column (x in axes widths), `text` ending at the top of the axes."""
    ax.annotate(text, (x, 1.0), xycoords='axes fraction', xytext=(0, 2), textcoords='offset points', va='bottom',
                fontsize=fs(6), color=GREY, linespacing=1.0, **kw)


def outside_figures(results, out):
    """Supplementary Figure 3 (si_fig_outside_data): Finch trained on ArchipelaGO, on each method's training split and on
    that split screened against the method's test proteins, and the method retrained by us, on all test proteins and on
    those without a relative; Supplementary Figure 4 (si_fig_outside_methods): Finch on the method's data against the
    method retrained by us and the authors' model, all test proteins."""
    folder = results / 'outside_masked'
    rows, groups = row_positions()
    unit = 0.215  # inches per row
    axes_h = unit * (rows[0][3] + HEADER + 0.45 - rows[-1][3] + 0.6)
    left_in, bottom_in, top_in = 1.08, 0.42, 0.5  # inches: the label column, the x axis, the key and column heads

    panel_w, column, gap = 2.2, 0.37, 0.16
    height = bottom_in + axes_h + top_in
    fig = plt.figure(figsize=(DOUBLE, height))
    panels = [fig.add_axes([(left_in + i * (panel_w + 2 * column + gap)) / DOUBLE, bottom_in / height, panel_w / DOUBLE,
                            axes_h / height]) for i in range(2)]
    for ax, subset, title in ((panels[0], 'all', 'All test proteins'),
                              (panels[1], 'no_relative', "No relative in the method's training split")):
        dot_panel(fig, ax, folder, rows, groups, subset, OUTSIDE_X, names=ax is panels[0])
        ax.set_title(title, fontsize=fs(6.5), fontweight='normal', loc='left', pad=3)
        own, screened_x = 1 + column / panel_w, 1 + 2 * column / panel_w  # right ends of the number columns, in axes widths
        column_head(ax, own, 'own\ndata', ha='right')
        column_head(ax, screened_x, 'screened\ndata', ha='right')
        ax.annotate('ArchipelaGO minus', ((own + screened_x) / 2 - column / panel_w / 2 + 0.01, 1.0), xycoords='axes fraction',
                    xytext=(0, 16), textcoords='offset points', ha='center', va='bottom', fontsize=fs(6), color=GREY)
        width = dot_width(fig, ax)
        for (method, label, test_set), version, _, y in rows:
            cell, diff = outside_cell(folder, method, test_set, subset, version, 'micro_AUPRC')
            v = {m: float(cell[m]['micro_AUPRC_mean']) for m in ('method', 'finch_A', 'finch_A_screened', 'finch_B') if m in cell}
            shift = dodge(v, width)
            ax.plot([v['finch_A'], v['finch_B']], [y + shift['finch_A'], y + shift['finch_B']], color=COLOURS['archipelago'],
                    lw=1.4, alpha=0.35, solid_capstyle='butt', zorder=2)
            ax.plot(v['method'], y + shift['method'], **method_style(method), zorder=3)
            ax.plot(v['finch_A'], y + shift['finch_A'], **FINCH_OWN, zorder=4)
            if 'finch_A_screened' in v:
                ax.plot(v['finch_A_screened'], y + shift['finch_A_screened'], **FINCH_SCREENED, zorder=4)
            ax.plot(v['finch_B'], y + shift['finch_B'], **FINCH_ARCHIPELAGO, zorder=5)
            screened = diff.get('finch_B - finch_A_screened')
            for x, text in ((own, signed(diff['finch_B - finch_A'])), (screened_x, signed(screened) if screened else '\u2013')):
                ax.text(x, y, text, transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=fs(6.5))
    key(fig, [("Finch on the method's data", FINCH_OWN), ('Finch on the screened data', FINCH_SCREENED),
              ('Finch on ArchipelaGO', FINCH_ARCHIPELAGO),
              ('Method on its data, retrained by us', [method_style(m) for m in ('deepgose', 'transfew', 'dpfunc')])])
    save(fig, out, 'si_fig_outside_data')

    width4, panel4 = 5.1, 3.3
    fig = plt.figure(figsize=(width4, height))
    ax = fig.add_axes([left_in / width4, bottom_in / height, panel4 / width4, axes_h / height])
    dot_panel(fig, ax, folder, rows, groups, 'all', (0.40, 0.80))
    column_head(ax, 1.02, 'Finch minus\nmethod', ha='left')
    width = dot_width(fig, ax)
    for (method, label, test_set), version, _, y in rows:
        cell, diff = outside_cell(folder, method, test_set, 'all', version, 'micro_AUPRC')
        v = {m: float(cell[m]['micro_AUPRC_mean']) for m in ('method', 'released', 'finch_A')}
        shift = dodge({m: v[m] for m in ('method', 'released')}, width)  # as Fig. 4b: Finch stays on the row
        ax.plot([v['method'], v['finch_A']], [y, y], color=COLOURS['archipelago'], lw=1.4, alpha=0.35, solid_capstyle='butt',
                zorder=2)
        ax.plot(v['method'], y + shift['method'], **method_style(method), zorder=3)
        ax.plot(v['released'], y + shift['released'], **method_style(method, hollow=True), zorder=3)
        ax.plot(v['finch_A'], y, **FINCH, zorder=4)
        ax.text(v['finch_A'], y + 0.2, r3(v['finch_A']), ha='center', va='bottom', fontsize=fs(6.5), zorder=5)
        ax.text(1.02, y, signed(diff['finch_A - method']), transform=ax.get_yaxis_transform(), ha='left', va='center',
                fontsize=fs(6.5))
    key(fig, [("Finch on the method's data", FINCH), ('Method, retrained by us', [method_style(m) for m in ('deepgose', 'transfew', 'dpfunc')]),
              ("Method, authors' model", [method_style(m, hollow=True) for m in ('deepgose', 'transfew', 'dpfunc')])])
    save(fig, out, 'si_fig_outside_methods')


def draw(results, out):
    out = out / 'si'
    out.mkdir(parents=True, exist_ok=True)
    tables(results, out)
    recipe_tables(results, out)
    cost_table(results, out)
    shared(results, out)
    outside_tables(results, out)
    vocabulary_table(results, out)
    outside_figures(results, out)
