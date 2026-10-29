"""Figure 3: choosing the computational proteins.

a, Finch's test micro-AUPRC gain over HQ alone (seed means) for the five ways of adding 150,000 proteins, the matrix of
ingredients each uses, and the paired step of each ingredient. Reads <root>/results/selection/{metrics.tsv,contrasts.tsv}.
b, how far the added proteins lie from HQ: the smoothed density of each added protein's distance to its nearest HQ protein
in ESM2-15B embedding space, for three of panel a's training sets (random computational proteins, random cluster
representatives, ArchipelaGO), and the median distance from an HQ protein to its nearest other HQ protein as the
reference (dashed: HQ's own density is a spike at zero). Reads <root>/results/distances/{densities.tsv,summary.tsv}
(retrain.scoring.distances).
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import Affine2D, ScaledTranslation, blended_transform_factory

from figures.style import DOUBLE, GREY, INK, contrasts, fs, half_up, read, save, seed_means

HQ, PICKED = '#4985d0', '#7565db'  # the HQ blue and ArchipelaGO violet of the other figures
H = 4.0  # figure height, inches
# panel a: training set, ingredients used (one per cluster, rare function, coverage ranking)
ARMS = [('hq_random', (0, 0, 0)), ('pool_random', (1, 0, 0)), ('filter_random', (1, 1, 0)),
        ('coverage_unfiltered', (1, 0, 1)), ('archipelago', (1, 1, 1))]
# row label and the paired step along random -> one per cluster -> + rare function -> + coverage (= ArchipelaGO)
STEPS = [('One per\ncluster', 'pool_random', 'hq_random'), ('Rare\nfunction', 'filter_random', 'pool_random'),
         ('Coverage\nranking', 'archipelago', 'filter_random')]
BAR, DOT_ON, DOT_OFF, STRIPE = '#8f939a', '#3b3e44', '#dcdde0', '#f3f3f1'
TOP, MS = 0.062, 6  # bar axis top; matrix dot size
# panel b: three of panel a's training sets, in its order, with line colour and key label (the two sets that differ only
# by the rare-function filter lie on the curves of the sets without it)
LINES = [('hq_random', '#c3c6cb', 'Random'), ('pool_random', '#7d8188', 'One per cluster'),
         ('archipelago', PICKED, 'ArchipelaGO')]


def at(ax, x, y, dx=0.0, dy=0.0):
    """A transform for drawing at (x, y) in axes-x and data-y coordinates, shifted by (dx, dy) inches."""
    base = blended_transform_factory(ax.transAxes, ax.transData)
    return (Affine2D() + ax.figure.dpi_scale_trans + ScaledTranslation(x, y, base)
            + ScaledTranslation(dx, dy, ax.figure.dpi_scale_trans))


def upset(b, m, gains, steps):
    """The bars of the gain over HQ alone (axes b) over the ingredient matrix (axes m, sharing b's x)."""
    for i, v in enumerate(gains):
        b.bar(i, v, 0.62, color=BAR, lw=0, zorder=2)
        b.text(i, v + TOP * 0.015, half_up(v, sign=True), ha='center', va='bottom', fontsize=fs(7.5), zorder=4)
    b.set_ylim(0, TOP)
    b.set_xlim(-0.6, len(ARMS) - 0.4)
    b.set_yticks([0, 0.02, 0.04, 0.06])
    b.set_ylabel('Micro-AUPRC gain\nover HQ alone')
    b.tick_params(axis='x', length=0, labelbottom=False)
    for row, (name, step) in enumerate(steps):  # rows top to bottom = the three ingredients
        y = len(steps) - 1 - row
        if row % 2 == 0:
            m.axhspan(y - 0.5, y + 0.5, color=STRIPE, lw=0, zorder=0)
        m.text(0, 0, name, transform=at(m, 0, y, dx=-0.05), ha='right', va='center', fontsize=fs(7.5), color=INK,
               linespacing=1.05, clip_on=False)
        m.text(0, 0, half_up(step, sign=True), transform=at(m, 1, y, dx=0.06), ha='left', va='center', fontsize=fs(7),
               color=GREY, clip_on=False)
    m.text(0, 0, 'Step', transform=at(m, 1, len(steps) - 0.5, dx=0.06, dy=0.03), ha='left', va='bottom', fontsize=fs(7),
           color=GREY, clip_on=False)
    for i, (_, used) in enumerate(ARMS):
        ys = [len(steps) - 1 - r for r, u in enumerate(used) if u]
        if len(ys) > 1:
            m.plot([i, i], [min(ys), max(ys)], color=DOT_ON, lw=1.6, zorder=2, solid_capstyle='round')
        for r, u in enumerate(used):
            m.plot(i, len(steps) - 1 - r, 'o', ms=MS, color=DOT_ON if u else DOT_OFF, mec='none', zorder=3, clip_on=False)
    m.text(len(ARMS) - 1, -0.62, 'ArchipelaGO', ha='center', va='top', fontsize=fs(7.5), color=INK)
    m.set_ylim(-0.5, len(steps) - 0.5)
    m.set_yticks([])
    m.tick_params(axis='x', length=0, labelbottom=False)
    for side in ('left', 'bottom'):
        m.spines[side].set_visible(False)


def densities(ax, results):
    density = {}
    for r in read(results / 'densities.tsv'):
        density.setdefault(r['condition'], []).append((float(r['distance']), float(r['density'])))
    hq_median = {r['condition']: float(r['median']) for r in read(results / 'summary.tsv')}['hq']
    ax.axvline(hq_median, color=HQ, lw=1.2, ls=(0, (3, 1.6)), label='HQ to nearest HQ (median)', zorder=2)
    for k, (name, color, text) in enumerate(LINES):
        x, y = np.array(density[name]).T
        ax.fill_between(x, y, color=color, alpha=0.22, lw=0, zorder=3 + k)
        ax.plot(x, y, color=color, lw=1.6 if name == 'archipelago' else 1.3, label=text, zorder=3 + k,
                solid_capstyle='butt')
    ax.spines['bottom'].set_zorder(10)  # the axis line over the filled densities
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, None)
    ax.set_xticks(np.arange(0, 0.51, 0.1))
    ax.set_xlabel('Distance to the nearest HQ protein')
    ax.set_ylabel('Density')
    handles, labels = ax.get_legend_handles_labels()
    order = list(range(1, len(handles))) + [0]  # the training sets in panel a's order, then the HQ reference
    ax.legend([handles[i] for i in order], [labels[i] for i in order], loc='upper right', frameon=False, fontsize=fs(6.5),
              handlelength=1.6, borderaxespad=0.1, labelspacing=0.35)


def letter(fig, x, text):
    fig.text(x / fig.get_figwidth(), 1 - 0.02 / H, text, ha='left', va='top', fontsize=fs(8.5), fontweight='bold')


def axes(fig, left, bottom, width, height, **kw):
    w = fig.get_figwidth()
    return fig.add_axes([left / w, bottom / H, width / w, height / H], **kw)


def panel_upset(fig, results, left, right):
    """Panel a between `left` (its letter) and `right` (its bars' right edge; the steps are printed beyond it)."""
    metric = seed_means(results)
    gains = [float(metric['finch', s, 'auprc']['value']) - float(metric['finch', 'hq', 'auprc']['value']) for s, _ in ARMS]
    differences = contrasts(results)
    steps = [(name, float(differences['auprc', 'finch', a, 'finch', b]['difference'])) for name, a, b in STEPS]
    bars_left = left + 0.62
    matrix_bottom, matrix_height = 0.2, 0.95
    bars_bottom = matrix_bottom + matrix_height + 0.16
    b = axes(fig, bars_left, bars_bottom, right - bars_left, H - 0.22 - 0.02 - bars_bottom)
    m = axes(fig, bars_left, matrix_bottom, right - bars_left, matrix_height, sharex=b)
    upset(b, m, gains, steps)


def draw(results, out):
    fig = plt.figure(figsize=(DOUBLE, H))
    letter(fig, 0, 'a')
    panel_upset(fig, results / 'selection', 0, 3.02)
    letter(fig, 3.49, 'b')
    densities(axes(fig, 3.96, 0.42, DOUBLE - 3.96 - 0.1, H - 0.42 - 0.24), results / 'distances')
    save(fig, out, 'fig3_choosing')
