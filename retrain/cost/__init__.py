"""Training cost of the published methods with their authors' original training code (Supplementary Note 9, the "original
implementation" columns of its table).

The paper's runs of TransFew and DPFunc used ports of their authors' training loops, and DeepGO-SE's a faster
reimplementation of its training step that computes the same objective (retrain.deepgose.train). Their cost with the
authors' own code was measured on one H100 per method with short timing jobs, and multiplied out to the epochs the paper's
runs trained:

    python -m retrain.cost.deepgose [--march AXIOMS]   the official DeepGO-SE step (deepgo2's FastTensorDataLoader, batch 37,
                                                       Adam, BCE plus every EL loss), ms per batch, on all ArchipelaGO proteins
    python -m retrain.cost.dpfunc                      DPFunc's released data path (graphs in RAM, GraphDataLoader) against
                                                       ours, on a random subset of ArchipelaGO, and its per-epoch validation
    python -m retrain.cost.transfew off|released       TransFew's released training loop, without and with its debugging
                                                       and logging switches, on a random subset of ArchipelaGO
    python -m retrain.cost.table                       the arithmetic of the table, from the measured times and the recorded
                                                       times of the paper's runs

Each harness prints its measurements (ms per batch, seconds per validation pass); table.py holds the values the paper's
timing jobs printed and the run times the paper's fits recorded, and reproduces the table. The harnesses' timed loops are the
code that ran; their data loading uses this repository's modules (the same inputs). A rerun measures the speed of the
hardware it runs on, so it agrees with the paper's times only on the same hardware and software.
"""
