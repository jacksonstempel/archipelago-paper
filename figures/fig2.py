"""Figure 2: where the benchmark and ArchipelaGO's proteins come from, as a Sankey diagram drawn to scale from the counts
recorded when the datasets were built (figures/sankey_counts.json).

Proteins flow down two trunks, HQ (a) and the computational candidates (b); bands leave where proteins are removed, and
the screen against the benchmark leaves both trunks at the same height. Band widths are to scale: column a and
ArchipelaGO at A_SCALE, column b (the rare-function filter included) ZOOM times narrower; the fan marks the change for the
selected proteins.

Geometry: every band is a closed polygon that reaches TUCK under the node bar or shape it meets and overlaps the trunk it
leaves by SEAM, so shapes never merely touch (a shared edge leaves an antialiasing seam). Vertical bands keep their
horizontal width along the S, the Sankey convention. Every curve is curvature-continuous (G2) and meets straight runs with
zero curvature: the S blends with a quintic, a turn eases its direction with the same quintic, and a band's edges are
offsets of its centreline. Colours change only at a node bar or smoothly: a removal and the stub of unselected candidates
fade to white; the selected band turns from the candidates' grey to violet, and the fan fades from violet towards the bar.

matplotlib draws a fade as an image clipped to the band, which its PDF stores as a raster. With --tikz, the PDF is
instead written through TikZ (figures/tikz.py, needs `tectonic` on PATH), which turns every fade into a true vector
gradient, as in the paper.
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon

from figures import tikz
from figures.drawing import C, GAP, box, canvas, flow, header
from figures.style import GREY, INK, fs, save, tint

COUNTS = Path(__file__).with_name('sankey_counts.json')
CANDIDATE = '#c4c8ce'
A_SCALE, ZOOM = 1.0 / 100_000, 100  # inches per protein; candidate-column reduction
A_BAND = tint(C['hq'], 0.6)
TUCK, SEAM = 0.012, 0.008  # inches a band reaches under what it meets; inches it overlaps the trunk it leaves
STRAIGHT = 0.06  # inches a band runs straight out of its trunk before it turns


def counts(path=COUNTS):
    """Protein counts; stops if any recorded balance does not hold."""
    record = json.loads(Path(path).read_text(encoding='utf-8'))
    n = {key: item['count'] for key, item in record['counts'].items()}
    for total, parts in record['balances']:
        if n[total] != sum(n[p] for p in parts):
            raise ValueError(f"{path}: {total} = {n[total]:,}, but {' + '.join(parts)} = {sum(n[p] for p in parts):,}")
    return n


def smooth(u):
    """Quintic smoothstep: 0 to 1 with zero first and second derivatives at both ends."""
    return u * u * u * (u * (6 * u - 15) + 10)


def shade(ax, poly, axis, a0, a1, c0, c1, z=1):
    """Fill `poly` with a straight gradient: c0 at coordinate a0, c1 at a1 along x (axis 0) or y (axis 1), constant beyond.
    Drawn as an image clipped to `poly`; the image keeps the fade, so tikz.save_pdf draws a true gradient."""
    poly = np.asarray(poly, float)
    lo, hi = poly.min(0), poly.max(0)
    u = np.clip((np.linspace(lo[axis], hi[axis], 512) - a0) / (a1 - a0), 0, 1)
    im = ax.imshow(u[None, :] if axis == 0 else u[:, None], cmap=LinearSegmentedColormap.from_list('shade', [c0, c1]),
                   vmin=0, vmax=1, extent=(lo[0], hi[0], lo[1], hi[1]), origin='lower', aspect='auto',
                   interpolation='bilinear', zorder=z)
    clip = Polygon(poly, closed=True, fc='none', ec='none', lw=0)
    ax.add_patch(clip)
    im.set_clip_path(clip)
    im.fade = (poly, axis, a0, a1, c0, c1)


def band(ax, top, bottom, color, z=1):
    """A vertical S-band from the x-interval `top` = (xl, xr, y) down to `bottom` = (xl, xr, y): each edge blends its x
    with the quintic while y runs linearly. Reaches TUCK beyond both ends."""
    (l0, r0, y0), (l1, r1, y1) = top, bottom
    t = np.linspace(0, 1, 320)
    edge = lambda x0, x1: np.column_stack([x0 + (x1 - x0) * smooth(t), y0 + (y1 - y0) * t])  # noqa: E731
    poly = np.vstack([[(l0, y0 + TUCK)], edge(l0, l1), [(l1, y1 - TUCK), (r1, y1 - TUCK)], edge(r0, r1)[::-1],
                      [(r0, y0 + TUCK)]])
    ax.add_patch(Polygon(poly, closed=True, fc=color, ec='none', lw=0, zorder=z))


def turn(theta0, theta1, R, m=320):
    """A quarter turn from direction theta0 to theta1 spanning R in x and in y, starting at the origin; the direction
    eases with the quintic (G2 with the straight runs it joins). Returns the points and the direction at each."""
    u = np.linspace(0, 1, m)
    th = theta0 + (theta1 - theta0) * smooth(u)
    d = np.column_stack([np.cos(th), np.sin(th)])
    p = np.vstack([[0.0, 0.0], np.cumsum((d[1:] + d[:-1]) / 2, axis=0) * (u[1] - u[0])])
    return p * (R / np.abs(p[-1]).max()), th


def offsets(points, theta, w):
    """The two edges of a band of width w whose centreline runs through `points` with direction `theta`."""
    n = np.column_stack([-np.sin(theta), np.cos(theta)])
    return points + n * w / 2, points - n * w / 2


def turn_size(w, R):
    """The turn span for a band of width w: at least 1.2 w, so the inner edge stays smooth."""
    return max(R, 1.2 * w)


def elbow(ax, xa, xb, y_top, y_c, x_end, color, fade_to=None, R=0.12):
    """A band of thickness xb - xa that leaves a trunk downward from y_top, turns through a G2 quarter turn spanning R and
    runs horizontally, centred on y_c, to x_end; with fade_to, it fades from `color` to fade_to along the horizontal run."""
    w, s = xb - xa, 1 if x_end > xb else -1
    xm, R = (xa + xb) / 2, turn_size(w, R)
    y_t = y_c + R  # where the turn begins
    if y_t > y_top:
        raise ValueError(f'elbow: no room for a turn of {R:.2f} in above y = {y_c:.2f}')
    p, th = turn(-np.pi / 2, 0.0 if s > 0 else -np.pi, R)
    p += (xm, y_t)
    right, left = offsets(p, th, w)
    outer, inner = (left, right) if s > 0 else (right, left)
    xo, xi = xm - s * w / 2, xm + s * w / 2
    poly = np.vstack([[(xo - s * SEAM, y_top + TUCK), (xo - s * SEAM, y_t)], outer, [(x_end, y_c - w / 2),
                      (x_end, y_c + w / 2)], inner[::-1], [(xi, y_t), (xi, y_top + TUCK)]])
    if fade_to is None:
        ax.add_patch(Polygon(poly, closed=True, fc=color, ec='none', lw=0, zorder=1))
    else:
        shade(ax, poly, 0, xm + s * R, x_end, color, fade_to, 1)


def zband(ax, x0, y0, x1, y1, y_h, width, color, R=0.12, z=2, fade_from=None):
    """A thin band from (x0, y0) down, along y_h and down again to (x1, y1), with G2 quarter turns spanning R; its
    thickness is `width` everywhere. With fade_from = (colour, xa, xb) it is that colour up to x = xa along its horizontal
    run and fades to `color` by x = xb."""
    s = 1 if x1 > x0 else -1
    t1, th1 = turn(-np.pi / 2, 0.0 if s > 0 else -np.pi, R)
    t2, th2 = turn(0.0 if s > 0 else -np.pi, -np.pi / 2, R)
    points = np.vstack([[(x0, y0 + TUCK)], t1 + (x0, y_h + R), t2 + (x1 - s * R, y_h), [(x1, y1 - TUCK)]])
    theta = np.concatenate([[-np.pi / 2], th1, th2, [-np.pi / 2]])
    a, b = offsets(points, theta, width)
    poly = np.vstack([a, b[::-1]])
    if fade_from is None:
        ax.add_patch(Polygon(poly, closed=True, fc=color, ec='none', lw=0, zorder=z))
    else:
        shade(ax, poly, 0, fade_from[1], fade_from[2], fade_from[0], color, z)


def strip(ax, x0, x1, top, bottom, color, z=1):
    """A filled rectangle drawn as a six-point polygon: PDF rasterisers snap the edges of a plain rectangle fill to whole
    pixels without antialiasing, which makes a node bar stand out from the antialiased band edges it sits on."""
    xm = (x0 + x1) / 2
    ax.add_patch(Polygon([(x0, bottom), (xm, bottom), (x1, bottom), (x1, top), (xm, top), (x0, top)], closed=True, fc=color,
                         ec='none', lw=0, zorder=z))


def removal(ax, xa, xb, y0, x_end, color, main, detail):
    """A band leaving a trunk's edge [xa, xb] at y0, turning outward to x_end and fading out; the count above it, the
    reason below."""
    half = (xb - xa) / 2
    y = y0 - STRAIGHT - turn_size(2 * half, 0.12)  # the horizontal run's centre
    elbow(ax, xa, xb, y0, y, x_end, color, fade_to='white')
    ha = 'right' if x_end > xb else 'left'
    ax.text(x_end, y + half + 0.03, main, ha=ha, va='bottom', fontsize=fs(7), color=INK)
    ax.text(x_end, y - half - 0.03, detail, ha=ha, va='top', fontsize=fs(6.3), color=GREY, linespacing=1.25)


def node(ax, x0, x1, y, color, name, detail):
    """A Sankey node: a short bar across the trunk at y, named on its left."""
    strip(ax, x0, x1, y + 0.018, y - 0.018, color, z=3)
    ax.text(x0 - 0.08, y + 0.012, name, ha='right', va='bottom', fontsize=fs(7.4), color=INK)
    ax.text(x0 - 0.08, y - 0.012, detail, ha='right', va='top', fontsize=fs(6.4), color=GREY)


def figure(n):
    # The 664 reconstructed GPCRdb variants are drawn as never having entered: they come off the specialist sources and
    # the reference instead of leaving as a removal of their own.
    n = dict(n, specialist=n['specialist'] - n['gpcrdb_variants'], reference=n['reference'] - n['gpcrdb_variants'])
    wa = {k: n[k] * A_SCALE for k in n}  # widths in column a and ArchipelaGO
    wb = {k: n[k] * A_SCALE / ZOOM for k in n}  # widths in column b

    W, HGT = 7.2, 4.35
    fig, ax = canvas(W, HGT)
    Y_SRC, Y_REF = HGT - 0.62, HGT - 1.12  # source nodes; tops of both trunks (the benchmark leaves under the Reference bar)
    Y_REL = Y_REF - 0.55  # relatives of benchmark proteins leave both trunks
    Y_NODE = Y_REL - 0.78  # HQ and candidates; the rare-function filter leaves under the Candidates bar
    Y_KEPT = Y_NODE - 0.56  # candidates that pass it; the selected proteins leave under its bar
    Y_END = Y_KEPT - 0.40  # the stub of unselected candidates ends; the fan's apex
    Y_STUB = Y_END - 0.22  # the stub of unselected candidates fades out by here
    Y_BAR, BAR = Y_END - 0.44, 0.27  # top and height of the ArchipelaGO bar
    header(ax, 0.03, HGT - 0.14, 'a', 'HQ')
    header(ax, 4.5, HGT - 0.14, 'b', 'Computational candidates')

    # a: sources (each protein counted once, under the first that applies) merge into the reference; neighbouring bands
    # overlap by SEAM under the Reference bar so no seam shows where they meet.
    left, right = 1.25, 1.25 + wa['reference']
    x, sources = left, []
    for key, name, shift in (('goa_experimental', 'GOA\nexperimental', 0.0), ('specialist', 'Specialist\ndatabases', 0.45),
                             ('goa_tas_ic', 'GOA\nTAS or IC', 0.9)):
        sources.append((key, name, x, x + wa[key], shift))
        x += wa[key]
    for i, (key, name, xl, xr, shift) in enumerate(sources):
        strip(ax, xl + shift, xr + shift, Y_SRC + 0.035, Y_SRC, C['hq'], z=3)
        band(ax, (xl + shift, xr + shift, Y_SRC),
             (xl - (SEAM if i else 0), xr + (SEAM if i < len(sources) - 1 else 0), Y_REF), A_BAND)
        ax.text((xl + xr) / 2 + shift, Y_SRC + 0.07, f'{n[key]:,}', ha='center', va='bottom', fontsize=fs(6.6), color=GREY)
        ax.text((xl + xr) / 2 + shift, Y_SRC + 0.2, name, ha='center', va='bottom', fontsize=fs(6.8), color=INK,
                linespacing=1.2)
    node(ax, left, right, Y_REF, C['hq'], 'Collected', f"{n['reference']:,}")
    # The trunk is one polygon: it steps in on the right under the Reference bar, where the benchmark leaves into its box,
    # and on the left at Y_REL, where the relatives leave; then straight on into ArchipelaGO, at the same scale.
    right1, left1 = right - wa['dev'] - wa['test'], left + wa['reference_related']
    ax.add_patch(Polygon([(left, Y_REF + TUCK), (right, Y_REF + TUCK), (right, Y_REF), (right1, Y_REF), (right1, Y_BAR - TUCK),
                          (left1, Y_BAR - TUCK), (left1, Y_REL), (left, Y_REL)], closed=True, fc=A_BAND, ec='none', lw=0,
                         zorder=1))
    BOX = (3.3, Y_REF - 0.55, 1.22, 0.56)  # centre x, centre y, width, height
    elbow(ax, right1, right, Y_REF, BOX[1], BOX[0] - BOX[2] / 2 + TUCK, A_BAND, R=0.2)
    box(ax, *BOX, 'Benchmark', f"{n['dev']:,} development\n{n['test']:,} test", 'white', GREY)
    removal(ax, left, left1, Y_REL, 0.03, A_BAND, f"{n['reference_related']:,} removed", 'related to the benchmark')
    node(ax, left1, right1, Y_NODE, C['hq'], 'HQ', f"{n['hq']:,}")

    # b: one random member of each sequence group enters the trunk, which steps in at each removal and, under the Filtered
    # candidates bar, on the left where the selected proteins leave; the rest ends as a stub.
    cl, cr = 5.0, 5.0 + wb['groups']
    cx = (cl + cr) / 2
    box(ax, cx, HGT - 0.54, 2.0, 0.44, 'TrEMBL, computational annotations only',
        f"{n['computational_pool']:,} of {n['trembl']:,} sequences")
    flow(ax, (cx, HGT - 0.76), (cx, Y_REF + 0.018 + GAP), tangent='v')
    ax.text(cx + 0.08, (HGT - 0.76 + Y_REF) / 2, 'one per 50% identity group', ha='left', va='center', fontsize=fs(6.3),
            color=GREY)
    node(ax, cl, cr, Y_REF, GREY, 'Representatives', f"{n['groups']:,}")
    # Not drawn: the 102 representatives with labels outside the vocabulary (too few to show; the caption says so), the
    # group members left unused, and the count of the proteins not selected (the stub says only that).
    cr1, cr2, xs = cr - wb['comp_related'], cr - wb['comp_related'] - wb['filter_removed'], cl + wb['selected']
    ax.add_patch(Polygon([(cl, Y_REF + TUCK), (cr, Y_REF + TUCK), (cr, Y_REL), (cr1, Y_REL), (cr1, Y_NODE), (cr2, Y_NODE),
                          (cr2, Y_KEPT), (cl, Y_KEPT)], closed=True, fc=CANDIDATE, ec='none', lw=0, zorder=1))
    # The stub of unselected candidates, then fading out: exactly the trunk minus the selected band. It reaches SEAM under
    # the band only while the band runs straight down beside it (to Y_TURN), so it never shows where the band turns away.
    Y_H = Y_KEPT - 0.2  # the selected band's horizontal run
    Y_TURN = Y_H + 0.12
    shade(ax, [(xs - SEAM, Y_KEPT + TUCK), (cr2, Y_KEPT + TUCK), (cr2, Y_STUB), (xs, Y_STUB), (xs, Y_TURN),
               (xs - SEAM, Y_TURN)], 1, Y_KEPT - 0.1, Y_STUB, CANDIDATE, 'white')
    removal(ax, cr1, cr, Y_REL, W - 0.1, CANDIDATE, f"{n['comp_related']:,} removed", 'related to the benchmark')
    removal(ax, cr2, cr1, Y_NODE, W - 0.1, CANDIDATE, f"{n['filter_removed']:,} removed", 'rare-function filter')
    node(ax, cl, cr1, Y_NODE, GREY, 'Candidates', f"{n['candidates']:,}")
    node(ax, cl, cr2, Y_KEPT, GREY, 'Filtered candidates', f"{n['filter_kept']:,}")
    ax.text((xs + cr2) / 2, Y_STUB - 0.04, 'not selected', ha='center', va='top', fontsize=fs(6.3), color=GREY)

    # ArchipelaGO, at column a's scale: the HQ proteins arrive straight. The selected proteins leave under the Filtered
    # candidates bar as a thin band that turns from the candidates' grey to violet, runs left and drops into the fan, which
    # enlarges them ZOOM times and fades towards the bar.
    v0, v1 = left1 + wa['hq'], left1 + wa['hq'] + wa['selected']
    apex = (v0 + v1) / 2
    xg = (cl + xs) / 2 - 0.12  # grey through the turn, then 0.9 in of fade
    zband(ax, (cl + xs) / 2, Y_KEPT, apex, Y_END, Y_H, wb['selected'], C['selected'], fade_from=(CANDIDATE, xg, xg - 0.9))
    t0, t1 = apex - wb['selected'] / 2, apex + wb['selected'] / 2
    # The fan meets the Selected bar exactly at its corners, then reaches TUCK straight down under it.
    fan = [(t0, Y_END), (t1, Y_END), (v1, Y_BAR), (v1, Y_BAR - TUCK), (v0, Y_BAR - TUCK), (v0, Y_BAR)]
    shade(ax, fan, 1, Y_END, Y_BAR, C['selected'], tint(C['selected'], 0.8))
    ax.text(apex, Y_BAR + 0.09, f'×{ZOOM}', ha='center', va='bottom', fontsize=fs(6.3), color=C['selected'], zorder=3)
    strip(ax, v0 - SEAM, v1, Y_BAR, Y_BAR - BAR, C['selected'], z=3)  # under the HQ part's edge: no seam between the two
    strip(ax, left1, v0, Y_BAR, Y_BAR - BAR, C['hq'], z=3)
    for x0, w, name, key in ((left1, wa['hq'], 'HQ', 'hq'), (v0, wa['selected'], 'Selected', 'selected')):
        ax.text(x0 + w / 2, Y_BAR - BAR / 2, f'{name}\n{n[key]:,}', ha='center', va='center', fontsize=fs(6.6),
                color='white', linespacing=1.15, zorder=4)
    ax.text(v1 + 0.14, Y_BAR - BAR / 2 + 0.01, 'ArchipelaGO', ha='left', va='bottom', fontsize=fs(8.5), fontweight='bold',
            color=INK)
    ax.text(v1 + 0.14, Y_BAR - BAR / 2 - 0.02, f"{n['archipelago']:,} proteins, {n['vocabulary']:,} GO terms",
            ha='left', va='top', fontsize=fs(6.6), color=GREY)
    if n['reference'] - n['dev'] - n['test'] - n['reference_related'] != n['hq']:
        raise ValueError('Sankey column a does not balance')
    return fig


def draw(out, use_tikz=False):
    """fig2_dataset.{png,pdf,svg} in `out`; with use_tikz the PDF is written through TikZ (true gradients)."""
    fig = figure(counts())
    if not use_tikz:
        save(fig, out, 'fig2_dataset')
        return
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / 'fig2_dataset.png')
    fig.savefig(out / 'fig2_dataset.svg', metadata={'Date': None})
    tikz.save_pdf(fig, out / 'fig2_dataset.pdf')
    plt.close(fig)
