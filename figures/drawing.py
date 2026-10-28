"""Drawing helpers of the schematic figures (Fig. 2 and Extended Data Fig. 1): connectors with arrowheads, boxes, panel
headings, a score strip. Coordinates are inches on an axes that spans the whole figure."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path as Curve

from figures.style import GREY, INK, fs

MONO = ['Consolas', 'DejaVu Sans Mono', 'Courier New']
C = {'hq': '#4985d0', 'selected': '#7565db'}  # HQ blue, ArchipelaGO violet
LINE, FILL = '#9da2a9', '#f4f5f6'
GAP = 0.035  # space between an arrowhead's tip and the thing it points at


def canvas(w, h):
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, w)
    ax.set_ylim(0, h)
    ax.axis('off')
    return fig, ax


def header(ax, x, y, letter, name):
    """Panel heading: bold letter, then the name."""
    renderer = ax.figure.canvas.get_renderer()
    t = ax.text(x, y, letter, ha='left', va='bottom', fontsize=fs(8.5), fontweight='bold', color=INK)
    right = ax.transData.inverted().transform((t.get_window_extent(renderer).x1, 0))[0]
    ax.text(right + 0.07, y, name, ha='left', va='bottom', fontsize=fs(8.5), color=INK)


def head(ax, tip, direction, color, lw=0.9):
    """A solid arrowhead with its tip at `tip`, pointing along `direction`; no shaft of its own."""
    d = np.asarray(direction, float)
    d /= np.linalg.norm(d)
    ax.add_patch(FancyArrowPatch(tuple(np.asarray(tip) - d * 0.045), tuple(tip), arrowstyle='-|>', mutation_scale=8, lw=lw,
                                 color=color, shrinkA=0, shrinkB=0, zorder=4))


def cubic(points, t):
    p0, p1, p2, p3 = points
    return (((1 - t) ** 3)[:, None] * p0 + (3 * (1 - t) ** 2 * t)[:, None] * p1 + (3 * (1 - t) * t ** 2)[:, None] * p2
            + (t ** 3)[:, None] * p3)


def first_part(points, t):
    """Control points of the cubic from its start to parameter t (de Casteljau)."""
    p0, p1, p2, p3 = points
    a, b, c = p0 + (p1 - p0) * t, p1 + (p2 - p1) * t, p2 + (p3 - p2) * t
    d, e = a + (b - a) * t, b + (c - b) * t
    return [p0, a, d, d + (e - d) * t]


def controls(p0, p1, tangent, k):
    """The two inner control points of a cubic from p0 to p1. Tangents: 'h' (horizontal at both ends), 'v' (vertical at
    both ends), 'vh' (leave vertically, arrive horizontally), 'hv' (the reverse)."""
    (x0, y0), (x1, y1) = p0, p1
    if tangent == 'h':
        return (x0 + (x1 - x0) * k, y0), (x1 - (x1 - x0) * k, y1)
    if tangent == 'vh':
        return (x0, y0 + (y1 - y0) * k), (x1 - (x1 - x0) * k, y1)
    if tangent == 'hv':
        return (x0 + (x1 - x0) * k, y0), (x1, y1 - (y1 - y0) * k)
    return (x0, y0 + (y1 - y0) * k), (x1, y1 - (y1 - y0) * k)


def flow(ax, p0, p1, color=LINE, lw=0.9, ls='-', arrow=True, tangent='h', k=0.5, z=1):
    """Smooth cubic connector from p0 to p1, with an arrowhead at p1; tangents as in controls()."""
    c1, c2 = controls(p0, p1, tangent, k)
    points = [np.asarray(q, float) for q in (p0, c1, c2, p1)]
    if arrow:
        t = np.linspace(0, 1, 2000)
        curve = cubic(points, t)
        distance = np.linalg.norm(curve - points[3], axis=1)
        cut = np.flatnonzero(distance >= 0.03)[-1]
        approach = curve[np.flatnonzero(distance >= 0.06)[-1]]  # the head's axis follows the line's last stretch
        points = first_part(points, t[cut])
    path = Curve([tuple(q) for q in points], [Curve.MOVETO, Curve.CURVE4, Curve.CURVE4, Curve.CURVE4])
    ax.add_patch(FancyArrowPatch(path=path, arrowstyle='-', lw=lw, color=color, linestyle=ls, capstyle='butt', zorder=z))
    if arrow:
        head(ax, np.asarray(p1, float), np.asarray(p1, float) - approach, color, lw)


def box(ax, x, y, w, h, name, detail='', face=FILL, edge=LINE, size=7.2, lw=0.8):
    """Rounded box centred on (x, y); the name and its one- or two-line detail are centred as a block."""
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0,rounding_size=0.07', fc=face, ec=edge,
                                lw=lw, zorder=2))
    if detail:
        shift = 0.052 * detail.count('\n')
        ax.text(x, y + 0.008 + shift, name, ha='center', va='bottom', fontsize=fs(size), color=INK, zorder=3)
        ax.text(x, y - 0.022 + shift, detail, ha='center', va='top', fontsize=fs(size - 1.0), color=GREY, linespacing=1.2,
                zorder=3)
    else:
        ax.text(x, y, name, ha='center', va='center', fontsize=fs(size), color=INK, zorder=3)


def scores(ax, x, y, w, h, color, seed):
    """A strip of bars standing for scores of many GO terms, a few of them high (schematic)."""
    values = np.random.default_rng(seed).beta(0.6, 6, 26) * 0.42
    values[[5, 13, 21]] = [0.94, 0.62, 0.47]
    bw = w / len(values)
    for i, v in enumerate(values):
        ax.add_patch(plt.Rectangle((x + i * bw + bw * 0.14, y), bw * 0.72, v * h, fc=color if v > 0.43 else '#c9ccd1',
                                   ec='none', zorder=2))
    ax.plot([x, x + w], [y, y], color=GREY, lw=0.5)
