"""Finch on any training set: train its readouts per seed (train), then their weighted mean (combine).

Fits go to runs/finch/<training set>/seed<k>/: the readouts in <readout>/ (readout.pt, raw dev_scores.npy and
test_scores.npy, history.json) and Finch itself (raw dev_scores.npy, test_scores.npy, record.json), which the scorer reads.
The models of the recipe ablation (Fig. 4a) combine readouts the same way into runs/<model>/archipelago/seed<k>/.
Recipes, weights and seeds come from configs/finch.json.
"""
import json
from pathlib import Path

CONFIG = Path(__file__).resolve().parents[2] / 'configs' / 'finch.json'


def config():
    return json.loads(CONFIG.read_text(encoding='utf-8'))
