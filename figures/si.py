"""Tables and figures of the Supplementary Information: LaTeX table bodies (rows only; the SI typesets the headers) and two
figures, written into <out>/si/. Reads, under <root>/results/:

  main/{metrics.tsv,contrasts.tsv}                         every model on every training set (retrain.scoring.score main)
  dataset_sizes.tsv                                        retrain.scoring.sizes
  experimental_test_labels/{metrics.tsv,contrasts.tsv,counts.tsv}   retrain.scoring.test_labels
  explabels_training/{metrics.tsv,contrasts.tsv}           retrain.scoring.score explabels_training
  no_hit/{subsets.tsv,differences.tsv}                     retrain.scoring.no_hit
  pr_curves/curves.tsv                                     retrain.scoring.pr_curves
  selection/{metrics.tsv,contrasts.tsv}                    retrain.scoring.score selection
  ladder/{metrics.tsv,contrasts.tsv}                       retrain.scoring.score ladder

Supplementary Table 1 (the annotation sources) is typed in the SI and the label-review prompt is released with the
reviewed labels; neither is generated here. Nor is the SI's comparison of the two losses on the development set, which
comes from the learning-rate and dropout search that chose the recipe (not part of this repository).
"""
import matplotlib.pyplot as plt
import numpy as np

from figures.style import COLOURS, DOUBLE, INK, NAMES, PUBLISHED, contrasts, fs, half_up, read, save, seed_means
from retrain.finch import config as finch_config

SETS = [('original', 'Original'), ('goa', 'GOA'), ('hq', 'HQ'), ('hq_random', 'HQ + random'), ('archipelago', 'ArchipelaGO')]
SETNAME = dict(SETS)
DASH = '--'


def r3(value, sign=False, places=3):
    """Rounded half up, as in the main text, with a LaTeX minus sign."""
    return half_up(value, sign, places).replace('-', '$-$')


def ci(low, high, places=3):
    return f'[{r3(low, places=places)}, {r3(high, places=places)}]'


def write(out, name, lines):
    with open(out / name, 'w', encoding='utf-8', newline='\n') as handle:
        handle.write('\n'.join(lines) + '\n')


def tables(results, out):
    """Every table body of the SI."""
    main = {split: seed_means(results / 'main', split) for split in ('test', 'development')}
    m = lambda method, ts, metric, split='test': main[split].get((method, ts, metric))  # noqa: E731
    contrast = contrasts(results / 'main')

    # Supplementary Table 2: composition of the training sets and the benchmark
    sizes = {r['dataset']: r for r in read(results / 'dataset_sizes.tsv')}
    lines = []
    for key, name in (('original_transfew', 'TransFew original data'), ('original_deepgose', 'DeepGO-SE original data'),
                      ('original_dpfunc', 'DPFunc original data'), ('goa', 'GOA'), ('hq', 'HQ'),
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
    def differences(table, method, a, method_b, b, metrics=(('auprc', 3), ('fmax', 3), ('smin', 2))):
        cells = []
        for k, places in metrics:
            r = table[k, method, a, method_b, b]
            cells.append(f"{r3(r['difference'], True, places)} {ci(r['ci_low'], r['ci_high'], places)}")
        return cells

    lines = []
    for method in PUBLISHED:
        for a, b in (('hq', 'original'), ('archipelago', 'hq'), ('archipelago', 'original')):
            lines.append(f'{NAMES[method]} & {SETNAME[a]} $-$ {SETNAME[b]} & '
                         + ' & '.join(differences(contrast, method, a, method, b)) + r' \\')
        lines.append(r'\addlinespace')
    for a, b in (('hq', 'goa'), ('hq_random', 'hq'), ('archipelago', 'hq'), ('archipelago', 'hq_random')):
        lines.append(f'Finch & {SETNAME[a]} $-$ {SETNAME[b]} & '
                     + ' & '.join(differences(contrast, 'finch', a, 'finch', b)) + r' \\')
    lines.append(r'\addlinespace')
    for ts in ('hq', 'archipelago'):
        for method in PUBLISHED:
            lines.append(f'Finch $-$ {NAMES[method]} & on {SETNAME[ts]} & '
                         + ' & '.join(differences(contrast, 'finch', ts, method, ts)) + r' \\')
    write(out, 'table_contrasts.tex', lines)
    lines = []
    for method in PUBLISHED:  # the GOA steps of the published methods
        for a, b in (('goa', 'original'), ('hq', 'goa')):
            lines.append(f'{NAMES[method]} & {SETNAME[a]} $-$ {SETNAME[b]} & '
                         + ' & '.join(differences(contrast, method, a, method, b)) + r' \\')
    write(out, 'table_contrasts_goa.tex', lines)

    # development against test
    lines = []
    for method in ['finch', 'transfew', 'deepgose', 'dpfunc', 'embedding_knn', 'diamond_knn', 'naive']:
        first = True
        for ts, name in SETS:
            t, d = m(method, ts, 'auprc'), m(method, ts, 'auprc', 'development')
            if t is None or d is None:
                continue
            lines.append(f"{NAMES[method] if first else ''} & {name} & {r3(d['value'])} & {r3(t['value'])} & "
                         f"{r3(float(t['value']) - float(d['value']), True)} \\\\")
            first = False
    write(out, 'table_devtest.tex', lines)

    # TransFew's checkpoint rule: last epoch (as published), lowest development BCE, highest development AUPRC
    lines = []
    variants = ('transfew', 'transfew_bcebest', 'transfew_auprcbest')
    for ts, name in SETS:
        if m('transfew', ts, 'auprc') is None:
            continue
        values = [m(v, ts, 'auprc') for v in variants]
        devs = [m(v, ts, 'auprc', 'development') for v in variants]
        lines.append(f'{name} & ' + ' & '.join(f"{r3(v['value'])} {ci(v['ci_low'], v['ci_high'])}" for v in values)
                     + ' & ' + ' & '.join(r3(v['value']) for v in devs) + r' \\')
    write(out, 'table_transfew.tex', lines)

    # experimental-only test labels
    folder = results / 'experimental_test_labels'
    counts = {r['variant']: r for r in read(folder / 'counts.tsv')}
    experimental = seed_means(folder)
    experimental_contrast = contrasts(folder)
    value = {'official': lambda method, ts: float(m(method, ts, 'auprc')['value']),
             'exp': lambda method, ts: float(experimental[method, ts, 'auprc']['value'])}
    models = ([(method, ts) for method in PUBLISHED for ts in ('original', 'goa', 'hq', 'archipelago')]
              + [('finch', ts) for ts in ('goa', 'hq', 'hq_random', 'archipelago')]
              + [('naive', ts) for ts in ('goa', 'hq', 'archipelago')])
    lines = []
    for method in ['finch', 'transfew', 'deepgose', 'dpfunc', 'naive']:
        first = True
        for ts, name in SETS:
            if (method, ts) not in models:
                continue
            lines.append(f"{NAMES[method] if first else ''} & {name} & {r3(value['official'](method, ts))} & "
                         f"{r3(value['exp'](method, ts))} \\\\")
            first = False
    write(out, 'table_explabels.tex', lines)
    steps = [(method, a, b) for method in PUBLISHED
             for a, b in (('original', 'goa'), ('goa', 'hq'), ('hq', 'archipelago'), ('original', 'archipelago'))] + \
            [('finch', 'goa', 'hq'), ('finch', 'hq', 'hq_random'), ('finch', 'hq', 'archipelago'),
             ('finch', 'hq_random', 'archipelago')]
    lines = []
    for method, a, b in steps:
        cells = []
        for labels, table in (('official', contrast), ('exp', experimental_contrast)):
            r = table['auprc', method, b, method, a]
            cells.append(f"{r3(value[labels](method, b) - value[labels](method, a), True)} {ci(r['ci_low'], r['ci_high'])}")
        lines.append(f'{NAMES[method]} & {SETNAME[b]} $-$ {SETNAME[a]} & ' + ' & '.join(cells) + r' \\')
        if b == 'archipelago' and a == 'original':
            lines.append(r'\addlinespace')
    write(out, 'table_explabels_diff.tex', lines)
    with open(out / 'explabels_counts.tex', 'w', encoding='utf-8', newline='\n') as handle:
        for variant, key in (('all', 'All'), ('exp', 'Exp')):
            r = counts[variant]
            handle.write(f"\\newcommand{{\\Count{key}Proteins}}{{{int(r['proteins']):,}}}\n"
                         f"\\newcommand{{\\Count{key}Pairs}}{{{int(r['propagated_pairs']):,}}}\n"
                         f"\\newcommand{{\\Count{key}Clusters}}{{{int(r['clusters']):,}}}\n"
                         f"\\newcommand{{\\Count{key}Terms}}{{{int(r['terms_with_positive']):,}}}\n")
        handle.write(f"\\newcommand{{\\CountExpDropped}}{{{int(counts['exp']['proteins_dropped']):,}}}\n")

    # HQ proteins trained with their experimental labels only
    folder = results / 'explabels_training'
    trained = seed_means(folder)
    trained_contrast = contrasts(folder)
    lines = []
    for method in ['finch'] + PUBLISHED:
        cells = []
        for ts in ('original', 'hq_explabels', 'hq', 'archipelago'):
            v = trained.get((method, ts, 'auprc'))
            cells.append(f"{r3(v['value'])} {ci(v['ci_low'], v['ci_high'])}" if v else DASH)
        lines.append(f'{NAMES[method]} & ' + ' & '.join(cells) + r' \\')
    write(out, 'table_ablate_hq.tex', lines)
    lines = []
    for method in ['finch'] + PUBLISHED:
        cells = []
        for a, b in (('hq_explabels', 'hq'), ('original', 'hq_explabels')):
            d = trained_contrast.get(('auprc', method, b, method, a))
            cells.append(f"{r3(d['difference'], True)} {ci(d['ci_low'], d['ci_high'])}" if d else DASH)
        lines.append(f'{NAMES[method]} & ' + ' & '.join(cells) + r' \\')
    write(out, 'table_ablate_hq_diff.tex', lines)

    # test proteins without training relatives
    folder = results / 'no_hit'
    subsets = [('none', 5000, 2734)] + [(f"{float(r['cutoff']):g}", int(r['proteins']), int(r['clusters']))
                                         for r in read(folder / 'subsets.tsv')]
    nohit = {(r['subset'], r['method'], r['base']): r for r in read(folder / 'differences.tsv')}
    cutname = {'none': 'all test proteins', '1e-10': r'$10^{-10}$', '1e-06': r'$10^{-6}$', '0.001': r'$10^{-3}$',
               '0.01': r'$10^{-2}$', '0.1': r'$10^{-1}$', '1': '1', '10': '10'}
    for name, columns in (('table_nohit.tex', [('finch', 'hq'), ('finch', 'hq_random')] +
                           [(method, 'original') for method in PUBLISHED]),
                          ('table_nohit_hq.tex', [(method, 'hq') for method in PUBLISHED])):
        lines = []
        for cut, n, k in subsets:
            cells = [f"{r3(nohit[cut, method, base]['difference'], True)} "
                     f"{ci(nohit[cut, method, base]['ci_low'], nohit[cut, method, base]['ci_high'])}"
                     for method, base in columns]
            lines.append(f'{cutname[cut]} & {n:,} & {k:,} & ' + ' & '.join(cells) + r' \\')
        write(out, name, lines)

    # selection rules and number of picks (Finch)
    folder = results / 'selection'
    selection = {split: seed_means(folder, split) for split in ('test', 'development')}
    selection_contrast = contrasts(folder)
    arms = [('hq', 'HQ (no proteins added)'), ('hq_random', '+ 150,000 random computational proteins (HQ + random)'),
            ('pool_random', '+ 150,000 random cluster representatives'),
            ('filter_random', '+ 150,000 random representatives passing the rare-function filter'),
            ('coverage_unfiltered', '+ the first 150,000 of the coverage ranking, without the filter'),
            ('archipelago', '+ the first 150,000 of the coverage ranking (ArchipelaGO)'),
            ('archipelago_300k', '+ the first 300,000 of the coverage ranking'),
            ('archipelago_600k', '+ the first 600,000 of the coverage ranking')]
    lines = []
    for key, name in arms:
        t = lambda k: selection['test']['finch', key, k]  # noqa: E731
        cells = [f"{r3(t('auprc')['value'])} {ci(t('auprc')['ci_low'], t('auprc')['ci_high'])}"]
        cells += [r3(t(k)['value']) for k in ('fmax', 'f1')]
        cells += [f"{float(t('smin')['value']):.2f}", r3(selection['development']['finch', key, 'auprc']['value'])]
        if key == 'archipelago_300k':
            lines.append(r'\midrule')
        lines.append(f'{name} & ' + ' & '.join(cells) + r' \\')
    write(out, 'table_arms.tex', lines)
    lines = []
    for (a, b), name in ((('pool_random', 'hq_random'), 'One per cluster (random representatives $-$ HQ + random)'),
                         (('filter_random', 'pool_random'),
                          'Rare-function filter (filtered $-$ unfiltered random representatives)'),
                         (('archipelago', 'filter_random'),
                          'Coverage ranking (ArchipelaGO $-$ filtered random representatives)'),
                         (('coverage_unfiltered', 'pool_random'),
                          'Coverage ranking without the filter ($-$ random representatives)'),
                         (('archipelago', 'coverage_unfiltered'),
                          'Filter alongside the ranking (ArchipelaGO $-$ ranking alone)'),
                         (('archipelago', 'hq_random'), 'All three (ArchipelaGO $-$ HQ + random)'),
                         (('archipelago_300k', 'archipelago'), '300,000 $-$ 150,000 picks'),
                         (('archipelago_600k', 'archipelago'), '600,000 $-$ 150,000 picks'),
                         (('archipelago_600k', 'archipelago_300k'), '600,000 $-$ 300,000 picks')):
        if a == 'archipelago_300k':
            lines.append(r'\midrule')
        cells = differences(selection_contrast, 'finch', a, 'finch', b,
                            (('auprc', 3), ('fmax', 3), ('f1', 3), ('smin', 2)))
        lines.append(f'{name} & ' + ' & '.join(cells) + r' \\')
    write(out, 'table_arms_diff.tex', lines)
    return subsets, nohit


def recipe_table(results, out):
    """The Finch recipe ablation (Fig. 4a), all on ArchipelaGO: test AUPRC, Fmax and Smin, development AUPRC, the step
    from the model it builds on and the difference from DeepGO-SE (paired)."""
    folder = results / 'ladder'
    ladder = {split: seed_means(folder, split) for split in ('test', 'development')}
    ladder_contrast = contrasts(folder)
    readouts = finch_config()['readouts']
    models = finch_config()['ladder']

    def setting(model):
        """The learning rate and dropout of a one-readout model, in LaTeX."""
        r = readouts[next(iter(models[model]))]
        mantissa, exponent = f"{r['lr']:.0e}".split('e')
        lr = r'$10^{%d}$' % int(exponent) if mantissa == '1' else r'$%s\times10^{%d}$' % (mantissa, int(exponent))
        return lr + ('' if r['dropout'] is None else f", {r['dropout']:g}")

    lines = []
    for model, name, base, where in (
            ('ladder_linear_esm2_3b', 'Linear readout on ESM2-3B', None, None),
            ('ladder_linear_esm2_15b', 'Linear readout on ESM2-15B', 'ladder_linear_esm2_3b', None),
            ('ladder_narrow', '+ hidden layer, 1,024 units, no dropout', 'ladder_linear_esm2_15b', None),
            ('ladder_hidden', '+ width 8,192 and dropout', 'ladder_narrow', None),
            ('ladder_asam', '+ ASAM', 'ladder_hidden', None),
            ('finch', '+ ProtT5 network, averaged (Finch)', 'ladder_asam', 'Methods'),
            ('ladder_two_esm2_15b', 'Control: + a second ESM2-15B network, averaged', 'ladder_asam', 'as step 4')):
        if model == 'ladder_two_esm2_15b':
            lines.append(r'\midrule')
        t = lambda k: ladder['test'][model, 'archipelago', k]  # noqa: E731
        step = DASH
        if base:
            d = ladder_contrast['auprc', model, 'archipelago', base, 'archipelago']
            step = f"{r3(d['difference'], True)} {ci(d['ci_low'], d['ci_high'])}"
        d = ladder_contrast['auprc', model, 'archipelago', 'deepgose', 'archipelago']
        development = ladder['development'][model, 'archipelago', 'auprc']['value']
        lines.append(f"{name} & {where or setting(model)} & {r3(t('auprc')['value'])} "
                     f"{ci(t('auprc')['ci_low'], t('auprc')['ci_high'])} & {r3(t('fmax')['value'])} & "
                     f"{float(t('smin')['value']):.2f} & {r3(development)} & "
                     f"{step} & {r3(d['difference'], True)} {ci(d['ci_low'], d['ci_high'])} \\\\")
    write(out, 'table_recipe.tex', lines)


def figures(results, out, subsets, nohit):
    """The gain against the no-hit cutoff, and the precision-recall curves of the published methods."""
    main = seed_means(results / 'main')
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 2.5), sharey=True)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.27, top=0.9, wspace=0.08)
    x = np.arange(len(subsets))
    ticks = ['All\n5,000'] + [f"{('$10^{%d}$' % round(np.log10(float(s)))) if float(s) < 0.5 else s}\n{n:,}"
                              for s, n, _ in subsets[1:]]
    shades = {'transfew': '#55595f', 'deepgose': '#8f939a', 'dpfunc': '#bfc2c7'}
    panels = [(axes[0], [(method, 'original', shades[method], NAMES[method]) for method in PUBLISHED],
               'Published methods: ArchipelaGO $-$ original data'),
              (axes[1], [('finch', 'hq', COLOURS['archipelago'], 'ArchipelaGO $-$ HQ'),
                         ('finch', 'hq_random', COLOURS['hq_random'], 'ArchipelaGO $-$ (HQ + random)')],
               'Finch: gain from the selected proteins')]
    for ax, series, title in panels:
        for method, base, color, name in series:
            d = np.array([[float(nohit[c, method, base][k]) for k in ('difference', 'ci_low', 'ci_high')]
                          for c, _, _ in subsets])
            ax.fill_between(x, d[:, 1], d[:, 2], color=color, alpha=0.12, lw=0)
            ax.plot(x, d[:, 0], color=color, lw=1.3, marker='o', ms=2.5, label=name)
        ax.axhline(0, color=INK, lw=0.6)
        ax.set_xticks(x, ticks, fontsize=fs(5.5))
        ax.set_title(title, fontsize=fs(7), fontweight='normal', loc='left')
        ax.legend(frameon=False, fontsize=fs(6), loc='upper right')
        ax.set_xlabel('No ultra-sensitive DIAMOND hit to any training protein at E $\\leq$ (proteins left)',
                      fontsize=fs(6))
    axes[0].set_ylabel('Gain in test micro-AUPRC')
    axes[0].set_ylim(-0.01, 0.2)
    save(fig, out, 'si_fig_nohit')

    pr = {}
    for r in read(results / 'pr_curves' / 'curves.tsv'):
        if r['seed'] == 'mean':
            pr.setdefault((r['method'], r['training_set']), []).append((float(r['recall']), float(r['precision'])))
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, 2.75), sharey=True)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.15, top=0.9, wspace=0.1)
    for ax, method in zip(axes, PUBLISHED):
        for ts, name in (('archipelago', 'ArchipelaGO'), ('hq', 'HQ'), ('original', 'Original data')):
            recall, precision = np.array(pr[method, ts]).T
            ax.plot(recall, precision, color=COLOURS[ts], lw=1.2, label=f"{name}  {r3(main[method, ts, 'auprc']['value'])} | "
                                                                       f"{r3(main[method, ts, 'fmax']['value'])}")
            f1 = 2 * recall * precision / np.maximum(recall + precision, 1e-12)
            k = int(np.argmax(f1))
            ax.scatter(recall[k], precision[k], s=12, color=COLOURS[ts], edgecolors=INK, linewidths=0.4, zorder=4)
        ax.set_title(NAMES[method], fontsize=fs(7), fontweight='bold', loc='left')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.02)
        ax.set_xticks(np.arange(0, 1.01, 0.2))
        ax.set_xlabel('Recall')
        ax.legend(loc='lower left', title='AUPRC | Fmax', frameon=False, fontsize=fs(5.5), title_fontsize=fs(5.5),
                  alignment='left', handlelength=1.4)
    axes[0].set_ylabel('Precision')
    save(fig, out, 'si_fig_pr')


def draw(results, out):
    out = out / 'si'
    out.mkdir(parents=True, exist_ok=True)
    subsets, nohit = tables(results, out)
    recipe_table(results, out)
    figures(results, out, subsets, nohit)
