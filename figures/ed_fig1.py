"""Extended Data Fig. 1: the Finch model, read left to right (schematic). The sequence goes to two frozen protein language
models, ESM2-15B (layer 48 of 48) and ProtT5 (layer 21 of 24, first 1,024 tokens); each one's residue vectors are averaged
into one vector, read by its own one-hidden-layer network (hidden width 8,192 and 4,096); the raw scores are averaged with
weights 0.6 and 0.4 and propagated up the GO graph to scores for the 8,739 molecular-function terms."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle

from figures.drawing import C, GAP, LINE, MONO, box, flow, scores
from figures.style import GREY, INK, fs, save, tint

VIOLET = C['selected']
W, H = 6.69, 1.98  # full text width (170 mm)
ROWS = {'ESM2-15B': 1.28, 'ProtT5': 0.44}  # row centres
MID = (ROWS['ESM2-15B'] + ROWS['ProtT5']) / 2
Y_HEAD = 1.68  # column headings (bottom)
SEQUENCE = 'MKTAYIAKQRQISF…'
# column x positions (inches): sequence, encoder box, pooled vector, readout network, weight, sum, GO graph, scores
X_SEQ, X_ENC, X_VEC, X_NET, X_W, X_SUM, X_GO, X_SC = 0.45, 1.62, 2.46, 3.12, 3.90, 4.38, 5.0, 6.1
ENC_W, ENC_H = 0.86, 0.46
VEC_W, VEC_H, VEC_N = 0.08, 0.52, 20
NET_W, NET_H = 0.72, 0.52
SC_W, SC_H = 1.0, 0.36


def heading(ax, x, text):
    ax.text(x, Y_HEAD, text, ha='center', va='bottom', fontsize=fs(7), color=INK, linespacing=1.15)


def vector(ax, x, y, seed):
    """The mean-pooled embedding: one column of cells shaded by value."""
    values = np.random.default_rng(seed).permutation(np.linspace(0, 1, VEC_N))  # the full range, light to dark
    h = VEC_H / VEC_N
    for i, v in enumerate(values):
        ax.add_patch(Rectangle((x - VEC_W / 2, y + VEC_H / 2 - (i + 1) * h), VEC_W, h, fc=tint(VIOLET, 0.94 - 0.94 * v),
                               ec='none', zorder=2))
    ax.add_patch(Rectangle((x - VEC_W / 2, y - VEC_H / 2), VEC_W, VEC_H, fill=False, ec=GREY, lw=0.5, zorder=3))


def network(ax, x, y, hidden):
    """Fully connected readout drawn left to right: input, one hidden layer (`hidden` nodes), output."""
    xs = np.linspace(x - NET_W / 2, x + NET_W / 2, 3)
    layers = [np.linspace(y - NET_H / 2 + 0.06, y + NET_H / 2 - 0.06, n) for n in (4, hidden, 4)]
    for xa, xb, ya, yb in zip(xs[:-1], xs[1:], layers[:-1], layers[1:]):
        for a in ya:
            for b in yb:
                ax.plot([xa, xb], [a, b], color=tint(VIOLET, 0.55), lw=0.3, zorder=1)
    for xx, ys in zip(xs, layers):
        ax.scatter(np.full(len(ys), xx), ys, s=8, color='white', edgecolor=VIOLET, lw=0.6, zorder=3)


def go_graph(ax, x, y):
    """A small GO graph, root on top: the leaf scored high passes its score up to its ancestors (violet)."""
    nodes = {'root': (0, 0.30), 'a': (-0.2, 0.05), 'b': (0.2, 0.05), 'l1': (-0.3, -0.22), 'l2': (-0.06, -0.22),
             'l3': (0.24, -0.22)}
    edges = [('a', 'root'), ('b', 'root'), ('l1', 'a'), ('l2', 'a'), ('l2', 'b'), ('l3', 'b')]
    hot = {'l2', 'a', 'b', 'root'}
    at = {k: np.array((x + dx, y + dy)) for k, (dx, dy) in nodes.items()}
    r = 0.045
    for child, parent in edges:
        on = child in hot and parent in hot
        p0, p1 = at[child], at[parent]
        d = (p1 - p0) / np.linalg.norm(p1 - p0)
        flow(ax, tuple(p0 + d * (r + 0.01)), tuple(p1 - d * (r + 0.012)), color=VIOLET if on else LINE,
             lw=0.8 if on else 0.6, tangent='v', k=0.0, z=1)
    for k, p in at.items():
        ax.add_patch(Circle(p, r, fc=VIOLET if k in hot else 'white', ec=VIOLET if k in hot else LINE, lw=0.7, zorder=3))


def figure():
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis('off')
    for x, text in ((X_SEQ, 'Protein\nsequence'), (X_ENC, 'Frozen language\nmodel'), (X_VEC, 'Mean over\nresidues'),
                    (X_NET, 'Readout\nnetwork'), (X_W, 'Weight'), (X_GO, 'Propagate\nup GO'),
                    (X_SC, 'Scores for 8,739\nGO terms')):
        heading(ax, x, text)
    seq = ax.text(0.06, MID, SEQUENCE, ha='left', va='center', fontsize=fs(6.5), family=MONO, color=INK)
    fig.canvas.draw()
    seq_right = ax.transData.inverted().transform((seq.get_window_extent().x1, 0))[0]
    details = {'ESM2-15B': 'layer 48 of 48', 'ProtT5': 'layer 21 of 24\nfirst 1,024 tokens'}
    for (name, y), weight, hidden, units in zip(ROWS.items(), (0.6, 0.4), (6, 4), ('8,192', '4,096')):
        fork = seq_right + 0.2  # the fan bends over a short stretch, then runs straight into the box
        flow(ax, (seq_right + 0.05, MID), (fork, y), tangent='h', k=0.5, arrow=False)
        flow(ax, (fork, y), (X_ENC - ENC_W / 2 - GAP, y))
        box(ax, X_ENC, y, ENC_W, ENC_H, name, details[name], size=7.0)
        flow(ax, (X_ENC + ENC_W / 2, y), (X_VEC - VEC_W / 2 - GAP, y))
        vector(ax, X_VEC, y, seed=len(name))
        flow(ax, (X_VEC + VEC_W / 2 + GAP, y), (X_NET - NET_W / 2 - 0.04 - GAP, y))
        network(ax, X_NET, y, hidden)
        ax.text(X_NET, y - NET_H / 2 - 0.05, f'{units} hidden', ha='center', va='top', fontsize=fs(6), color=GREY)
        flow(ax, (X_NET + NET_W / 2 + 0.04 + GAP, y), (X_W - 0.19 - GAP, y))
        ax.add_patch(FancyBboxPatch((X_W - 0.19, y - 0.1), 0.38, 0.2, boxstyle='round,pad=0,rounding_size=0.1',
                                    fc='white', ec=LINE, lw=0.7, zorder=2))
        ax.text(X_W, y, f'× {weight:.1f}', ha='center', va='center', fontsize=fs(6.5), color=INK, zorder=3)
        flow(ax, (X_W + 0.19 + GAP, y), (X_SUM, MID + np.sign(y - MID) * (0.075 + GAP)), tangent='hv', k=1.0)
    ax.add_patch(Circle((X_SUM, MID), 0.075, fc='white', ec=LINE, lw=0.7, zorder=2))
    ax.text(X_SUM, MID - 0.004, '+', ha='center', va='center', fontsize=fs(7), color=INK, zorder=3)
    go_graph(ax, X_GO, MID - 0.02)
    flow(ax, (X_SUM + 0.075 + GAP, MID), (X_GO - 0.36 - GAP, MID))
    flow(ax, (X_GO + 0.33 + GAP, MID), (X_SC - SC_W / 2 - GAP, MID))
    scores(ax, X_SC - SC_W / 2, MID - SC_H / 2, SC_W, SC_H, VIOLET, 3)
    return fig


def draw(out):
    save(figure(), out, 'ed_fig1_finch_model')
