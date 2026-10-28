# Third-party notices

This repository is under the MIT License (`LICENSE`). Parts of the retraining pipeline contain code copied or adapted
from the published methods it retrains. Their licenses permit redistribution, and their texts are kept next to the code.

| Code in this repository | Upstream | License | License text |
|---|---|---|---|
| `retrain/deepgose/model.py` (copied), `retrain/deepgose/axioms/makeEl.groovy` (copied), `retrain/deepgose/axioms/Normalizer.groovy` (adapted) | DeepGO-SE, [bio-ontology-research-group/deepgo2](https://github.com/bio-ontology-research-group/deepgo2), commit `6c77c985ef9453548dca629830c044b2dc1b4982` | BSD 3-Clause | `retrain/deepgose/LICENSE-deepgo2` |
| `retrain/transfew/model.py`, `retrain/transfew/label_gcn.py` (copied or adapted) | TransFew, [BioinfoMachineLearning/TransFew](https://github.com/BioinfoMachineLearning/TransFew), commit `5a4c9d7d726c95891deb7f7789622749a44bf62a` | MIT | `retrain/transfew/LICENSE-TransFew` |
| `retrain/dpfunc/model.py` (copied), `retrain/dpfunc/train.py`, `retrain/dpfunc/inputs.py` (adapted) | DPFunc, [CSUBioGroup/DPFunc](https://github.com/CSUBioGroup/DPFunc), commit `c5faf4df4e29b7d264ccb14d88254a136b54f4b2` | MIT | `retrain/dpfunc/LICENSE-DPFunc` |
| `retrain/dpfunc/embed.py` (`representation`, adapted from `esm.model.esm1.ProteinBertModel.forward`); `retrain/deepgose/features.py` (deepgo2's `extract_esm` batching on fair-esm) | fair-esm, [facebookresearch/esm](https://github.com/facebookresearch/esm) | MIT | `retrain/dpfunc/LICENSE-fair-esm` |

Each copied file names its upstream file and commit in its docstring or header.

## Software and data used but not included

These are installed or downloaded by the user; nothing of them is redistributed here.

| Item | Used for | License or terms |
|---|---|---|
| ESM2 and ESM-1b ([facebookresearch/esm](https://github.com/facebookresearch/esm), and the Hugging Face ports) | embeddings | MIT |
| ProtT5-XL-U50 ([Rostlab](https://huggingface.co/Rostlab/prot_t5_xl_half_uniref50-enc)) | Finch's second encoder | ProtTrans: Academic Free License 3.0 |
| ColabFold 1.6.3 ([sokrypton/ColabFold](https://github.com/sokrypton/ColabFold)) | structures for DPFunc | MIT; AlphaFold2 code Apache 2.0, AlphaFold2 parameters CC BY 4.0 |
| AlphaFold Protein Structure Database | structures for DPFunc | CC BY 4.0 |
| InterProScan 5.78-109.0 | DPFunc's domain inputs | Apache 2.0 (member databases have their own terms) |
| DIAMOND 2.2.4 | homology screen, DIAMOND-kNN | GPL-3.0 (called as an external program) |
| `pritamdeka/S-PubMedBert-MS-MARCO` | TransFew's label embeddings | see its model card |
| Gene Ontology, GO-Plus | labels, DeepGO-SE's axioms | CC BY 4.0 |
| The published methods' training data (DeepGO-SE, TransFew and DPFunc releases) | their original-data training sets | the authors' terms |

When you use these components, cite them as their authors ask.
