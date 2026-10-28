"""Shared figure style: Nature-family sizes (single column 89 mm, double column 183 mm), text at most 7 pt, bold lowercase
panel letters, the paper's colours, rounding as in the text, reading result tables and saving."""
import csv
from decimal import ROUND_HALF_UP, Decimal

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402

SINGLE, DOUBLE = 3.5, 7.2  # column widths in inches
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7, 'axes.titlesize': 8, 'axes.titleweight': 'bold', 'axes.labelsize': 7, 'xtick.labelsize': 6,
    'ytick.labelsize': 6, 'legend.fontsize': 6, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.axisbelow': True, 'savefig.dpi': 300, 'pdf.fonttype': 42, 'svg.fonttype': 'path', 'svg.hashsalt': 'fixed',
})
INK, GREY = '#1d1f22', '#5f636a'
# training sets: original data, HQ, HQ + random, ArchipelaGO
COLOURS = {'original': '#8f939a', 'hq': '#4985d0', 'hq_random': '#6fc3a0', 'archipelago': '#7565db'}
NAMES = {'transfew': 'TransFew', 'deepgose': 'DeepGO-SE', 'dpfunc': 'DPFunc', 'finch': 'Finch', 'naive': 'Naive',
         'embedding_knn': 'Embedding-kNN', 'diamond_knn': 'DIAMOND-kNN'}
PUBLISHED = ['transfew', 'deepgose', 'dpfunc']


def fs(points):
    """A font size, at most 7 pt."""
    return min(points, 7.0)


def tint(color, amount=0.55):
    """The colour mixed with white; amount 1 gives white."""
    r, g, b = to_rgb(color)
    return (r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount)


def half_up(value, sign=False, places=3):
    """Rounded half up (0.0195 -> 0.020), three decimals unless stated, as everywhere in the paper."""
    s = str(Decimal(repr(float(value))).quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP))
    if s.startswith('-') and float(s) == 0:
        s = s[1:]
    return '+' + s if sign and not s.startswith('-') else s


def read(path):
    with open(path, encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def seed_means(results, split='test'):
    """{(method, training set, metric): row} of the seed-mean rows of a results folder's metrics.tsv."""
    return {(r['method'], r['training_set'], r['metric']): r for r in read(results / 'metrics.tsv')
            if r['seed'] == 'mean' and r['split'] == split}


def contrasts(results):
    """{(metric, method, training set, minus method, minus training set): row} of a results folder's contrasts.tsv."""
    return {(r['metric'], r['method'], r['training_set'], r['minus_method'], r['minus_training_set']): r
            for r in read(results / 'contrasts.tsv')}


def pr_curves(results):
    """{(method, training set): [(recall, precision)]} of the `mean` curves (their area is the seed mean of AUPRC)."""
    curves = {}
    for r in read(results / 'pr_curves' / 'curves.tsv'):
        if r['seed'] == 'mean':
            curves.setdefault((r['method'], r['training_set']), []).append((float(r['recall']), float(r['precision'])))
    return curves


def save(fig, folder, name):
    """<name>.png, .pdf and .svg at the figure's exact size; no dates in the files, so an unchanged figure gives
    unchanged files."""
    folder.mkdir(parents=True, exist_ok=True)
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(folder / f'{name}.{ext}', metadata={'pdf': {'CreationDate': None}, 'svg': {'Date': None}}.get(ext))
    plt.close(fig)
