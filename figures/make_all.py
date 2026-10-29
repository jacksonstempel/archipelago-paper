"""Every figure and table of the paper from the result files under the data root.

    python -m figures.make_all [--data DIR] [--out DIR] [--tikz]

Reads <root>/results/ (written by retrain.scoring; see figures/__init__.py for the files each figure needs) and writes
into --out (default <root>/figures/): fig1_dataset, fig2_retraining, fig3_choosing, fig4_finch, ed_fig1_finch_model
(.png, .pdf, .svg), table1.tex, and si/ (the Supplementary Information's table bodies and figures). --tikz writes
Fig. 1's PDF through TikZ with true gradients, as in the paper (needs `tectonic` on PATH).
"""
import argparse
from pathlib import Path

from retrain.paths import add_argument, data_root
from figures import ed_fig1, fig1, fig3, main as main_figures, si


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argument(parser)
    parser.add_argument('--out', type=Path, help='output folder (default: <root>/figures)')
    parser.add_argument('--tikz', action='store_true', help="Fig. 1's PDF through TikZ and Tectonic")
    a = parser.parse_args()
    root = data_root(a.data)
    results, out = root / 'results', a.out or root / 'figures'
    out.mkdir(parents=True, exist_ok=True)
    fig1.draw(out, a.tikz)
    main_figures.figures(results, out)
    fig3.draw(results, out)
    ed_fig1.draw(out)
    si.draw(results, out)
    print('wrote', out)


if __name__ == '__main__':
    main()
