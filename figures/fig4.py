"""Figure 4: choosing the computational proteins.

a, the coverage ranking in four steps (a schematic, not data): in a two-dimensional stand-in for embedding space, each
candidate has a dashed line to its closest covered protein (HQ or an earlier pick); each tile shows the state after the
picks so far (violet dots) with the longest line, that of the next pick, in violet. The picks are not placed by hand: they
follow from the drawn positions by the rule itself (farthest_first).
b, Finch's test micro-AUPRC gain over HQ alone (seed means) for the five ways of adding 150,000 proteins, the matrix of
ingredients each uses, and the paired step of each ingredient. Reads <root>/results/selection/{metrics.tsv,contrasts.tsv}.
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import Affine2D, ScaledTranslation, blended_transform_factory

from figures.style import DOUBLE, GREY, INK, contrasts, fs, half_up, save, seed_means

HQ, PICKED, CAND = '#4985d0', '#7565db', '#a9adb3'  # the HQ blue and ArchipelaGO violet of the other figures
GAP_LINE = '#b9bcc2'  # a candidate's dashed line to its closest covered protein
DASH = (0, (2.2, 1.6))
# Panel a, in unit coordinates of the field, placed by hand: 14 HQ proteins around the edges and eight candidates between
# them, placed so that each pick's line is clearly the longest.
HQ_XY = np.array([(0.10, 0.90), (0.20, 0.94), (0.13, 0.80), (0.84, 0.95), (0.92, 0.88), (0.15, 0.05), (0.05, 0.20),
                  (0.93, 0.20), (0.80, 0.18), (0.05, 0.55), (0.10, 0.44), (0.86, 0.48), (0.48, 0.95), (0.55, 0.05)])
CAND_XY = np.array([(0.45, 0.50), (0.38, 0.36), (0.68, 0.73), (0.27, 0.30), (0.24, 0.84), (0.72, 0.30), (0.22, 0.62),
                    (0.78, 0.62)])
N_PICK = 4  # one tile per pick: each shows the state before it, with its line violet
H, A_LEFT, A_WIDTH = 4.0, 0.07, 3.68  # inches: figure height, panel a's left edge and width
# training set, ingredients used (one per cluster, rare function, coverage ranking)
ARMS = [('hq_random', (0, 0, 0)), ('pool_random', (1, 0, 0)), ('filter_random', (1, 1, 0)),
        ('coverage_unfiltered', (1, 0, 1)), ('archipelago', (1, 1, 1))]
# row label and the paired step along random -> one per cluster -> + rare function -> + coverage (= ArchipelaGO)
STEPS = [('One per\ncluster', 'pool_random', 'hq_random'), ('Rare\nfunction', 'filter_random', 'pool_random'),
         ('Coverage\nranking', 'archipelago', 'filter_random')]
BAR, DOT_ON, DOT_OFF, STRIPE = '#8f939a', '#3b3e44', '#dcdde0', '#f3f3f1'
TOP, MS = 0.062, 6  # bar axis top; matrix dot size
TILE, HGAP, VGAP, GRID_TOP = 1.68, 0.12, 0.24, 0.27  # tile side, gaps, space above the grid (inches)
S_HQ, S_CAND, PAD = 16, 14, 0.04  # dot areas (pt^2); margin of the field around the points
TITLES = ['Find farthest candidate', 'Pick and remeasure', 'Pick and remeasure', 'Repeat']
# each tile's labels: text, the point it names (an HQ protein, a candidate, or the n-th pick) and the direction of the
# leader from the point (degrees)
LABELS = [[('HQ protein', ('hq', 8), 180), ('candidate\nprotein', ('cand', 4), 0), ('farthest\nfrom HQ', ('pick', 0), 0)],
          [('picked', ('pick', 0), 0), ('new farthest', ('pick', 1), 180)],
          [('picked', ('pick', 1), 180), ('new farthest', ('pick', 2), 270)],
          []]
DOT_GAP, LEADER, TEXT_GAP = 1.6, 7.0, 1.6  # leader: gap from the dot's edge, length, gap to the text (points)


def farthest_first(candidates, covered, picks):
    """The coverage ranking: repeatedly pick the candidate farthest from everything covered (HQ and earlier picks)."""
    order = []
    for _ in range(picks):
        d = np.linalg.norm(candidates[:, None] - covered[None], axis=2)
        i = int(np.argmax(d.min(1)))  # earlier picks are covered, so their own distance is 0
        order.append(i)
        covered = np.vstack([covered, candidates[i]])
    return order


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


def letter(fig, x, text):
    fig.text(x / fig.get_figwidth(), 1 - 0.02 / H, text, ha='left', va='top', fontsize=fs(8.5), fontweight='bold')


def axes(fig, left, bottom, width, height, **kw):
    w = fig.get_figwidth()
    return fig.add_axes([left / w, bottom / H, width / w, height / H], **kw)


def label(ax, text, xy, angle, s):
    """A label for the dot of area s at xy: a straight leader along angle, stopping DOT_GAP short of the dot's edge, and
    the text beyond it."""
    u = np.array([np.cos(np.radians(angle)), np.sin(np.radians(angle))]).round(6)
    start = np.sqrt(s) / 2 + DOT_GAP  # the marker's diameter in points is sqrt(s)
    ax.annotate('', xy, xytext=tuple(u * (start + LEADER)), textcoords='offset points', xycoords='data',
                arrowprops=dict(arrowstyle='-', lw=0.5, color=GREY, shrinkA=0, shrinkB=start))
    ha = 'left' if u[0] > 0 else 'right' if u[0] < 0 else 'center'
    va = 'bottom' if u[1] > 0 else 'top' if u[1] < 0 else 'center'
    ax.annotate(text, xy, xytext=tuple(u * (start + LEADER + TEXT_GAP)), textcoords='offset points', ha=ha, va=va,
                fontsize=fs(6.5), color=INK, linespacing=1.1, zorder=5)


def tile(ax, hq, cand, picked, upcoming, labels=()):
    """The state after the picks so far: every other candidate has a dashed line to its closest covered protein; the line
    of the upcoming pick, the longest, is violet."""
    covered = np.vstack([hq] + [cand[i][None] for i in picked])
    for i in range(len(cand)):
        if i in picked:
            continue
        j = covered[np.argmin(np.linalg.norm(covered - cand[i], axis=1))]
        hot = i == upcoming
        ax.plot(*zip(cand[i], j), color=PICKED if hot else GAP_LINE, lw=1.4 if hot else 0.75, ls=DASH,
                dash_capstyle='round', zorder=2)
    rest = [i for i in range(len(cand)) if i not in picked]
    ax.scatter(*cand[rest].T, s=S_CAND, color=CAND, lw=0, zorder=3)
    ax.scatter(*hq.T, s=S_HQ, color=HQ, lw=0, zorder=3)
    if picked:
        ax.scatter(*cand[picked].T, s=S_HQ, color=PICKED, lw=0, zorder=3)
    for text, xy, angle, s in labels:
        label(ax, text, xy, angle, s)
    ax.set_xlim(-PAD, 1 + PAD)
    ax.set_ylim(-PAD, 1 + PAD)
    ax.axis('off')
    ax.add_patch(plt.Rectangle((-PAD, -PAD), 1 + 2 * PAD, 1 + 2 * PAD, facecolor='white', edgecolor='black', lw=0.5,
                               zorder=0, clip_on=False))


def panel_a(fig):
    hq, cand = HQ_XY, CAND_XY
    order = farthest_first(cand, hq, N_PICK)
    letter(fig, 0, 'a')
    for k in range(N_PICK):  # tile k: after k picks, with the line of pick k + 1 in violet
        col, row = k % 2, k // 2
        left = A_LEFT + col * (TILE + HGAP)
        bottom = H - GRID_TOP - TILE - row * (TILE + VGAP)
        ax = axes(fig, left, bottom, TILE, TILE)
        named = {'hq': lambda n: (hq[n], S_HQ), 'cand': lambda n: (cand[n], S_CAND),
                 'pick': lambda n: (cand[order[n]], S_HQ if n < k else S_CAND)}  # picked: violet, the size of HQ
        labels = [(text, named[kind](n)[0], angle, named[kind](n)[1]) for text, (kind, n), angle in LABELS[k]]
        tile(ax, hq, cand, order[:k], order[k], labels)
        ax.set_title(TITLES[k], loc='left', fontsize=fs(7), fontweight='normal', color=INK, pad=3)


def panel_b(fig, results):
    metric = seed_means(results)
    gains = [float(metric['finch', s, 'auprc']['value']) - float(metric['finch', 'hq', 'auprc']['value']) for s, _ in ARMS]
    differences = contrasts(results)
    steps = [(name, float(differences['auprc', 'finch', a, 'finch', b]['difference'])) for name, a, b in STEPS]
    left = A_WIDTH + 0.12
    letter(fig, left, 'b')
    bars_left, right = left + 0.62, DOUBLE - 0.42
    matrix_bottom, matrix_height = 0.2, 0.95
    bars_bottom = matrix_bottom + matrix_height + 0.16
    b = axes(fig, bars_left, bars_bottom, right - bars_left, H - 0.22 - 0.02 - bars_bottom)
    m = axes(fig, bars_left, matrix_bottom, right - bars_left, matrix_height, sharex=b)
    upset(b, m, gains, steps)


def draw(results, out):
    fig = plt.figure(figsize=(DOUBLE, H))
    panel_a(fig)
    panel_b(fig, results / 'selection')
    save(fig, out, 'fig4_choosing')
