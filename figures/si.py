"""Tables and figures of the Supplementary Information: LaTeX table bodies (rows only; the SI typesets the headers) and
one figure, written into <out>/si/. Reads, under <root>/results/:

  main/{metrics.tsv,contrasts.tsv}                         every model on every training set (retrain.scoring.score main)
  dataset_sizes.tsv                                        retrain.scoring.sizes
  experimental_test_labels/{metrics.tsv,contrasts.tsv,counts.tsv}   retrain.scoring.test_labels experimental
  all_goa_test_labels/{metrics.tsv,contrasts.tsv}          retrain.scoring.test_labels all_goa
  masked_test_labels/{metrics.tsv,contrasts.tsv}           retrain.scoring.masked
  explabels_training/{metrics.tsv,contrasts.tsv}           retrain.scoring.score explabels_training
  no_hit/{subsets.tsv,differences.tsv,E0.001/,E10/}        retrain.scoring.no_hit
  selection/{metrics.tsv,contrasts.tsv}                    retrain.scoring.score selection
  ladder/{metrics.tsv,contrasts.tsv}                       retrain.scoring.score ladder
  parameters.json                                          retrain.scoring.parameters
  training_cost.tsv                                        retrain.scoring.training_cost

Main-text Table 2 with intervals (table2_intervals.tex) is written beside Table 2 by figures/main.py.

Supplementary Table 1 (the annotation sources) is typed in the SI and the label-review prompt is released with the
reviewed labels; neither is generated here. Nor is the SI's comparison of the two losses on the development set, which
comes from the learning-rate and dropout search that chose the recipe (not part of this repository).
"""
import json

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
    def differences(table, method, a, method_b, b):
        cells = []
        for k, places in (('auprc', 3), ('fmax', 3), ('smin', 2)):
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

    experimental_contrast = label_sets(results, out)

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
    lines = []
    for cut, n, k in subsets:
        cells = [f"{r3(nohit[cut, method, base]['difference'], True)} "
                 f"{ci(nohit[cut, method, base]['ci_low'], nohit[cut, method, base]['ci_high'])}"
                 for method, base in [('finch', 'hq'), ('finch', 'hq_random')] + [(p, 'original') for p in PUBLISHED]]
        lines.append(f'{cutname[cut]} & {n:,} & {k:,} & ' + ' & '.join(cells) + r' \\')
    write(out, 'table_nohit.tex', lines)

    # Finch's lead over DeepGO-SE and TransFew, all trained on ArchipelaGO, on two of the no-hit subsets and with
    # experimental-only labels
    size = {cut: f'{n:,}' for cut, n, _ in subsets}
    lines = []
    for table, name, n in ((contrasts(folder / 'E0.001'), r'No hit at \Eval{10^{-3}}', size['0.001']),
                           (contrasts(folder / 'E10'), r'No hit at \Eval{10}', size['10']),
                           (experimental_contrast, 'Experimental-only labels', r'\CountExpProteins{}')):
        for method in ('deepgose', 'transfew'):
            cells = [paired(table[k, 'finch', 'archipelago', method, 'archipelago'], places)
                     for k, places in (('auprc', 3), ('fmax', 3), ('smin', 2))]
            lines.append((f'{name} & {n}' if method == 'deepgose' else ' & ') + f' & {NAMES[method]} & '
                         + ' & '.join(cells) + r' \\')
    write(out, 'table_lead.tex', lines)

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
        lines.append(f'{name} & ' + ' & '.join(differences(selection_contrast, 'finch', a, 'finch', b)) + r' \\')
    write(out, 'table_arms_diff.tex', lines)
    return subsets, nohit


def label_sets(results, out):
    """The gains between training sets with the reviewed (official), masked, experimental-only and all GOA test labels, in
    that column order: the difference of seed means in test micro-AUPRC with its paired 95% interval
    (table_explabels_diff.tex); and the counts of the experimental-only labels (explabels_counts.tex). Returns the
    experimental-only labels' contrasts."""
    folders = ('main', 'masked_test_labels', 'experimental_test_labels', 'all_goa_test_labels')
    value = {f: {(m, t): float(r['value']) for (m, t, k), r in seed_means(results / f).items() if k == 'auprc'}
             for f in folders}
    table = {f: contrasts(results / f) for f in folders}
    steps = [(method, a, b) for method in PUBLISHED
             for a, b in (('original', 'goa'), ('goa', 'hq'), ('hq', 'archipelago'), ('original', 'archipelago'))] + \
            [('finch', 'goa', 'hq'), ('finch', 'hq', 'hq_random'), ('finch', 'hq', 'archipelago'),
             ('finch', 'hq_random', 'archipelago')]
    lines = []
    for method, a, b in steps:
        cells = []
        for f in folders:
            r = table[f]['auprc', method, b, method, a]
            cells.append(f"{r3(value[f][method, b] - value[f][method, a], True)} {ci(r['ci_low'], r['ci_high'])}")
        lines.append(f'{NAMES[method]} & {SETNAME[b]} $-$ {SETNAME[a]} & ' + ' & '.join(cells) + r' \\')
        if b == 'archipelago' and a == 'original':
            lines.append(r'\addlinespace')
    write(out, 'table_explabels_diff.tex', lines)
    counts = {r['variant']: r for r in read(results / 'experimental_test_labels' / 'counts.tsv')}
    with open(out / 'explabels_counts.tex', 'w', encoding='utf-8', newline='\n') as handle:
        for variant, key in (('all', 'All'), ('exp', 'Exp')):
            r = counts[variant]
            handle.write(f"\\newcommand{{\\Count{key}Proteins}}{{{int(r['proteins']):,}}}\n"
                         f"\\newcommand{{\\Count{key}Pairs}}{{{int(r['propagated_pairs']):,}}}\n"
                         f"\\newcommand{{\\Count{key}Clusters}}{{{int(r['clusters']):,}}}\n"
                         f"\\newcommand{{\\Count{key}Terms}}{{{int(r['terms_with_positive']):,}}}\n")
        handle.write(f"\\newcommand{{\\CountExpDropped}}{{{int(counts['exp']['proteins_dropped']):,}}}\n")
    return table['experimental_test_labels']


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
    """The Finch recipe ablation (Fig. 4a) and further models, all trained on ArchipelaGO: test AUPRC, Fmax and Smin,
    development AUPRC, for the ablation the paired gain over the step it builds on, and trainable parameters in
    millions; then the paired comparisons cited in the main text."""
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
            ('ladder_two_esm2_15b', 'Control: + a second ESM2-15B network, averaged', 'ladder_asam', 'as step 4')):
        if model == 'ladder_two_esm2_15b':
            lines.append(r'\midrule')
        step = DASH
        if base:
            d = ladder_contrast['auprc', model, 'archipelago', base, 'archipelago']
            step = f"{r3(d['difference'], True)} {ci(d['ci_low'], d['ci_high'])}"
        lines.append(f"{name} & {where or setting(next(iter(models[model])))} & {scores(model)} & {step} & "
                     f"{millions(sum(parameters[model].values()))} \\\\")
    write(out, 'table_recipe.tex', lines)

    # further models; Finch with ProtT5 layer 24 shows the learning rate and dropout of that network
    lines = []
    for model, name, where in (('ladder_hidden_esm2_3b', 'Wide hidden layer on ESM2-3B', None),
                               ('ladder_hidden_2048', 'Hidden layer of 2,048 units on ESM2-15B', None),
                               ('finch_prott5_l24', 'Finch with ProtT5 layer 24', setting('prott5_l24')),
                               ('deepgose', 'DeepGO-SE', 'published'), ('transfew', 'TransFew', 'published'),
                               ('dpfunc', 'DPFunc', 'published')):
        parts = parameters[model].values()
        count = millions(sum(parts))
        if model == 'deepgose':  # one part per world of the ensemble
            lines.append(r'\midrule')
            count = f'{millions(sum(parts) / len(parts))} per model, {count} for six'
        lines.append(f"{name} & {where or setting(next(iter(models[model])))} & {scores(model)} & {count} \\\\")
    write(out, 'table_variants.tex', lines)

    # the paired comparisons cited in the main text
    lines = []
    for a, b, name in (('ladder_linear_esm2_15b', 'deepgose', 'Linear readout on ESM2-15B $-$ DeepGO-SE'),
                       ('ladder_narrow', 'deepgose', 'Hidden layer of 1,024 units on ESM2-15B $-$ DeepGO-SE'),
                       ('ladder_narrow', 'transfew', 'Hidden layer of 1,024 units on ESM2-15B $-$ TransFew'),
                       ('ladder_hidden', 'deepgose', 'Wide hidden layer on ESM2-15B $-$ DeepGO-SE'),
                       ('ladder_hidden_esm2_3b', 'deepgose', 'Wide hidden layer on ESM2-3B $-$ DeepGO-SE'),
                       ('deepgose', 'ladder_linear_esm2_3b', 'DeepGO-SE $-$ linear readout on ESM2-3B'),
                       ('ladder_hidden_esm2_3b', 'ladder_linear_esm2_3b',
                        'Wide hidden layer $-$ linear readout, both on ESM2-3B'),
                       ('ladder_hidden', 'transfew', 'Wide hidden layer on ESM2-15B $-$ TransFew'),
                       ('transfew', 'ladder_linear_esm2_15b', 'TransFew $-$ linear readout on ESM2-15B'),
                       ('ladder_hidden_2048', 'ladder_hidden', '2,048 units $-$ 8,192 units, both on ESM2-15B'),
                       ('finch', 'finch_prott5_l24', 'Finch $-$ Finch with ProtT5 layer 24')):
        cells = [paired(ladder_contrast[k, a, 'archipelago', b, 'archipelago'], places)
                 for k, places in (('auprc', 3), ('fmax', 3), ('smin', 2))]
        lines.append(f'{name} & ' + ' & '.join(cells) + r' \\')
    write(out, 'table_pairs.tex', lines)


def figure(out, subsets, nohit):
    """The gain against the no-hit cutoff."""
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


def draw(results, out):
    out = out / 'si'
    out.mkdir(parents=True, exist_ok=True)
    subsets, nohit = tables(results, out)
    recipe_tables(results, out)
    cost_table(results, out)
    figure(out, subsets, nohit)
