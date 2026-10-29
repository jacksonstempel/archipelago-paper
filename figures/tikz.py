"""Save a drawn matplotlib figure as a fully vector PDF through TikZ, so colour fades become true PDF gradients (matplotlib can
only fill a shape with one colour; shade() in figures/fig1.py fakes a fade with a clipped image, which this replaces).

The figure is drawn by matplotlib as usual; save_pdf() then reads the finished artists (outlines in display coordinates,
each text line's baseline as matplotlib placed it) and writes them as TikZ, which Tectonic compiles with Arial, the font
matplotlib used. Handles what the diagrams use: filled or stroked patches, text, and images carrying a `fade` attribute
(polygon, axis, a0, a1, c0, c1: c0 at coordinate a0, c1 at a1, constant beyond); anything else stops with an error.
"""
import io
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_pdf import PdfFile, RendererPdf
from matplotlib.colors import to_rgba
from matplotlib.image import AxesImage
from matplotlib.patches import FancyArrowPatch, Patch
from matplotlib.path import Path as MPath
from matplotlib.text import Text

SPECIAL = {'\\': r'\textbackslash{}', '{': r'\{', '}': r'\}', '%': r'\%', '&': r'\&', '#': r'\#', '_': r'\_', '$': r'\$',
           '~': r'\textasciitilde{}', '^': r'\textasciicircum{}'}


def colour(c):
    r, g, b, a = to_rgba(c)
    return f'{{rgb,1:red,{r:.4f};green,{g:.4f};blue,{b:.4f}}}', a


def path_tikz(path, scale):
    """A matplotlib path (display coordinates) as a TikZ path in inches."""
    out, start, last = [], None, None
    for verts, code in path.iter_segments(simplify=False, curves=True):
        p = np.asarray(verts, float).reshape(-1, 2) / scale
        xy = lambda q: f'({q[0]:.5f},{q[1]:.5f})'  # noqa: E731
        if code == MPath.MOVETO:
            out.append(xy(p[0])); start = last = p[0]
        elif code == MPath.LINETO:
            out.append('-- ' + xy(p[0])); last = p[0]
        elif code == MPath.CURVE3:  # quadratic: raise to a cubic
            c, end = p
            out.append(f'.. controls {xy(last + 2 / 3 * (c - last))} and {xy(end + 2 / 3 * (c - end))} .. {xy(end)}'); last = end
        elif code == MPath.CURVE4:
            out.append(f'.. controls {xy(p[0])} and {xy(p[1])} .. {xy(p[2])}'); last = p[2]
        elif code == MPath.CLOSEPOLY:
            out.append('-- cycle'); last = start
    return ' '.join(out)


def stroke_options(artist):
    edge, alpha = colour(artist.get_edgecolor())
    cap = {'butt': 'butt', 'round': 'round', 'projecting': 'rect'}[artist.get_capstyle()]
    return f'draw={edge}, draw opacity={alpha:.3f}, line width={artist.get_linewidth():.3f}bp, line cap={cap}, ' \
           f'line join={artist.get_joinstyle()}'


def patch_tikz(patch, scale):
    if isinstance(patch, FancyArrowPatch):
        paths, fillable = patch._get_path_in_displaycoord()
        face = patch.get_edgecolor()  # an arrow's head is filled in its line colour
    else:
        paths, fillable = [patch.get_transform().transform_path(patch.get_path())], [patch.get_fill()]
        face = patch.get_facecolor()
    lines = []
    stroked = patch.get_linewidth() > 0 and to_rgba(patch.get_edgecolor())[3] > 0
    assert patch.get_linestyle() in ('-', 'solid'), f'tikz: unsupported line style {patch.get_linestyle()!r}'
    for path, fill in zip(paths, fillable):
        options = []
        c, a = colour(face)
        if fill and a > 0:
            options.append(f'fill={c}, fill opacity={a:.3f}')
        if stroked:
            options.append(stroke_options(patch))
        if options:
            lines.append(rf'\path[{", ".join(options)}] {path_tikz(path, scale)};')
    return lines


def fade_tikz(image, scale, k):
    """The image's fade as a PDF axial shading of its polygon: \\pgfshadepath maps the middle 50 bp of a 100 bp shading onto
    the path's bounding box, so the colour at each position follows the fade exactly, constant beyond a0 and a1."""
    poly, axis, a0, a1, c0, c1 = image.fade
    data = image.axes.transData
    poly_px = data.transform(np.asarray(poly, float))
    lo, hi = np.asarray(poly, float).min(0)[axis], np.asarray(poly, float).max(0)[axis]
    c0, c1 = np.array(to_rgba(c0)[:3]), np.array(to_rgba(c1)[:3])
    at = lambda v: c0 + (c1 - c0) * np.clip((v - a0) / (a1 - a0), 0, 1)  # noqa: E731  the fade at a coordinate
    stops = sorted({lo, hi, *[a for a in (a0, a1) if lo < a < hi]})
    spec = '; '.join(f'rgb({25 + 50 * (v - lo) / (hi - lo):.4f}bp)=({",".join(f"{x:.4f}" for x in at(v))})'
                     for v in [stops[0]] + stops)  # the first stop twice: the shading's lower margin
    spec = spec.replace(f'rgb({25:.4f}bp)', 'rgb(0bp)', 1) + f'; rgb(100bp)=({",".join(f"{x:.4f}" for x in at(hi))})'
    kind = 'horizontal' if axis == 0 else 'vertical'
    path = MPath(np.vstack([poly_px, poly_px[:1]]), closed=True)  # a closed path's last vertex is only the close marker
    return [rf'\pgfdeclare{kind}shading{{fade{k}}}{{100bp}}{{{spec}}}',
            rf'\path[shading=fade{k}, shading angle=0] {path_tikz(path, scale)};']


def text_tikz(text, renderer, scale):
    """Each line of a text at the baseline matplotlib placed it on, left-aligned there (matplotlib has already applied the
    horizontal and vertical alignment)."""
    if not text.get_visible() or not text.get_text():
        return []
    assert text.get_rotation() == 0, 'tikz: rotated text is not supported'
    _, info, _ = text._get_layout(renderer)
    x0, y0 = text.get_transform().transform(text.get_unitless_position())
    c, a = colour(text.get_color())
    size = text.get_fontsize()
    weight = r'\bfseries' if text.get_fontweight() in ('bold', 700) else ''
    lines = []
    for line, _, dx, dy in info:
        s = ''.join(SPECIAL.get(ch, ch) for ch in line)
        lines.append(rf'\node[anchor=base west, inner sep=0, outer sep=0, text={c}, text opacity={a:.3f}, '
                     rf'font=\fontsize{{{size:.2f}bp}}{{{size * 1.2:.2f}bp}}\selectfont{weight}] '
                     rf'at ({(x0 + dx) / scale:.5f},{(y0 + dy) / scale:.5f}) {{{s}}};')
    return lines


def save_pdf(fig, path):
    """Write `fig` to `path` (PDF) through TikZ and Tectonic. The figure is first laid out as matplotlib lays out a PDF (72
    dpi, PDF renderer), not as on screen, where text heights are rounded to whole pixels."""
    w, h = fig.get_size_inches()
    dpi = fig.dpi
    fig.set_dpi(72)
    try:
        renderer = RendererPdf(PdfFile(io.BytesIO()), 72, h * 72, w * 72)
        renderer.file.newPage(w, h)
        fig.draw(renderer)
        body = tikz_body(fig, renderer, 72)
    finally:
        fig.set_dpi(dpi)
    doc = '\n'.join([r'\documentclass[border=0pt]{standalone}', r'\usepackage{fontspec}', r'\setmainfont{Arial}',
                     r'\usepackage{tikz}', r'\begin{document}', r'\begin{tikzpicture}[x=1in, y=1in]',
                     rf'\useasboundingbox (0,0) rectangle ({w:.5f},{h:.5f});', *body,
                     r'\end{tikzpicture}', r'\end{document}', ''])
    tectonic = shutil.which('tectonic')
    if not tectonic:
        raise SystemExit('tikz: the TikZ PDF needs Tectonic (https://tectonic-typesetting.github.io) on PATH')
    path = Path(path)
    with tempfile.TemporaryDirectory() as tmp:
        tex = Path(tmp) / (path.stem + '.tex')
        tex.write_text(doc, encoding='utf-8')
        run = subprocess.run([tectonic, '--outdir', tmp, str(tex)], capture_output=True, text=True,
                             env={**os.environ, 'SOURCE_DATE_EPOCH': '0'})  # no build time in the file: same figure, same bytes
        if run.returncode:
            raise RuntimeError(f'tikz: Tectonic failed on {path.name}:\n{run.stderr[-3000:]}')
        path.write_bytes((Path(tmp) / (path.stem + '.pdf')).read_bytes())


def tikz_body(fig, renderer, scale):
    """The figure's artists as TikZ commands, in matplotlib's drawing order; `scale` display units per inch."""
    body = []
    assert not fig.texts and not fig.patches and not fig.images, 'tikz: figure-level artists are not supported'
    k = 0
    for ax in fig.axes:
        for artist in sorted(ax.get_children(), key=lambda a: a.get_zorder()):  # stable, as matplotlib draws them
            if not artist.get_visible() or artist in ax.spines.values() or artist is ax.patch:
                continue
            if isinstance(artist, Text):
                body += text_tikz(artist, renderer, scale)
            elif isinstance(artist, Patch):
                body += patch_tikz(artist, scale)
            elif isinstance(artist, AxesImage):
                assert hasattr(artist, 'fade'), 'tikz: an image without a fade cannot be drawn as vector'
                body += fade_tikz(artist, scale, k)
                k += 1
            elif type(artist).__name__ in ('XAxis', 'YAxis'):
                assert not ax.axison, 'tikz: axes with ticks are not supported'
            else:
                raise TypeError(f'tikz: cannot draw {type(artist).__name__}')
    return body
