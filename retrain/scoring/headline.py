"""The paper's results on its main test labels, the annotations-masked version: <root>/results_masked/.

    python -m retrain.scoring.headline [--data DIR] [--device cuda|cpu]

The main test labels leave out of every metric the GOA assertions that the label review considered, accepted or not, and the
pairs they add through the ontology: neither positive nor negative (retrain.scoring.test_labels masked). Every result that
depends on the test labels is scored again from the same predictions with the scorer's --labels masked and written under
<root>/results_masked/, in the same folders with the same files as under <root>/results/ (figures/ reads them with
HEADLINE=masked):

  main, selection, ladder, explabels_training/    retrain.scoring.score RUN --labels masked
  no_hit/                                         retrain.scoring.no_hit --labels masked (the same subsets)
  term_strata/gains.tsv                           retrain.scoring.term_strata --labels masked
  macro/                                          retrain.scoring.macro --labels masked
  pr_curves/curves.tsv                            retrain.scoring.pr_curves --labels masked
  dataset_sizes.tsv                               retrain.scoring.sizes --labels masked (the test row)
  experimental_test_labels/counts.tsv             results/experimental_test_labels/counts.tsv with its `all` row giving the
                                                  masked labels' counts (from results/masked_test_labels/counts.tsv)

Needs the reviewed scoring of the four runs (scoring/RUN.json: python -m retrain.scoring.score RUN) and the label versions'
(python -m retrain.scoring.test_labels experimental and masked). A few minutes on one GPU.
"""
import argparse
import csv

from retrain.common import add_argument, data_root, log
from retrain.scoring import macro, no_hit, pr_curves, score, sizes, term_strata
from retrain.scoring.score import RUNS, results_dir
from retrain.scoring.test_labels import write_counts


def counts(root):
    """results_masked/experimental_test_labels/counts.tsv (see the module docstring)."""
    def read(path):
        with open(path, newline='') as handle:
            return list(csv.DictReader(handle, delimiter='\t'))

    reviewed = read(root / 'results' / 'experimental_test_labels' / 'counts.tsv')
    masked = next(r for r in read(root / 'results' / 'masked_test_labels' / 'counts.tsv') if r['variant'] == 'masked')
    rows = [dict(masked, variant='all') if r['variant'] == 'all' else r for r in reviewed]
    out = results_dir(root, 'masked') / 'experimental_test_labels' / 'counts.tsv'
    out.parent.mkdir(parents=True, exist_ok=True)
    write_counts(out, rows)
    log(f'wrote {out}')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--device', choices=['cuda', 'cpu'])
    add_argument(parser)
    a = parser.parse_args(argv)
    root = data_root(a.data)
    if results_dir(root, 'masked').exists():
        raise FileExistsError(f'{results_dir(root, "masked")} exists')
    options = ['--labels', 'masked', '--data', str(root)]
    device = ['--device', a.device] if a.device else []
    for run in RUNS:
        score.main([run, *options, *device])
    for module in (no_hit, term_strata, macro, pr_curves):
        module.main([*options, *device])
    sizes.main(options)
    counts(root)


if __name__ == '__main__':
    main()
