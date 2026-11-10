"""The arithmetic of the training-cost table (retrain.cost.table) gives the table of Supplementary Note 9."""
from retrain.cost.table import table


def test_the_training_cost_table_of_the_paper():
    text, final = table()
    rounded = {m: tuple(round(v, 2 if v < 1 else 1) for v in values) for m, values in final.items()}
    # (original implementation: one run on ArchipelaGO, all runs; as trained here: one run, all runs), H100-hours
    assert rounded == {'transfew': (0.73, 4.3, 0.63, 4.0), 'dpfunc': (3.7, 31.4, 3.4, 24.6),
                       'deepgose': (62.0, 143.1, 8.8, 20.7)}
    assert 'FINAL TABLE' in text
