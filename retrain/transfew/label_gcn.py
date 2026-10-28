"""TransFew's label graph and its pretrained GCN, over the output vocabulary on the paper's GO release.

Adapted from the official repository (MIT; LICENSE-TransFew), commit 5a4c9d7:
  - graph: DataGen/LabelEmbeddings. Node features are sentence embeddings of each term's text definition (the text
    between the quotes) from pritamdeka/S-PubMedBert-MS-MARCO (the paper says BioBERT; the code uses this model under
    the name "biobert"); edges are the direct child -> parent is_a and part_of relations among the terms.
  - pretraining: label_embedding.py with hparams.mf.label_ae: GAE(LabelEncoder(768, 1024)); NormalizeFeatures +
    RandomLinkSplit(num_val=0.2, num_test=0, is_undirected=False, split_labels=True, add_negative_train_samples=True);
    Adam lr 1e-4, weight decay 1e-2; 500 epochs; seed 42; the checkpoint with the lowest validation reconstruction
    loss (`<=`). No protein or label data are used.
The nodes are the vocabulary terms in vocabulary order (official: the terms with at least one training annotation),
so one graph and one GCN serve every training set.

    python -m retrain.transfew.label_gcn [--data DIR]       (one GPU, a few minutes)

Writes inputs/transfew/label_gcn/{graph.pt, gcn.pt}. Needs obonet, networkx, sentence-transformers, torch_geometric.
"""
import argparse
import json
import math
import random
from pathlib import Path

import numpy as np
import torch

from retrain import common
from retrain.transfew.model import LABEL_TEXT_DIM, label_gae

TEXT_MODEL = 'pritamdeka/S-PubMedBert-MS-MARCO'
TEXT_REVISION = '96786c7024f95c5aac7f2b9a18086c7b97b23036'


def folder(root):
    return Path(root) / 'inputs' / 'transfew' / 'label_gcn'


def seed_all(value):
    random.seed(value)
    np.random.seed(value)
    torch.manual_seed(value)
    torch.cuda.manual_seed_all(value)


def definitions_and_edges(obo, vocabulary):
    """Each term's definition text, and (child, parent) column pairs of the direct is_a and part_of edges."""
    import networkx as nx
    import obonet
    graph = obonet.read_obo(str(obo))
    graph.remove_edges_from([e for e in graph.edges(keys=True) if e[2] not in ('is_a', 'part_of')])
    if not nx.is_directed_acyclic_graph(graph):
        raise ValueError('the is_a/part_of graph has a cycle')
    missing = [t for t in vocabulary if t not in graph.nodes]
    if missing:
        raise ValueError(f'{len(missing)} vocabulary terms missing from the ontology, e.g. {missing[0]}')
    definitions = []
    for t in vocabulary:
        text = graph.nodes[t]['def']
        if text.count('"') != 2:   # the source's check: the definition is one quoted string
            raise ValueError(f'{t}: definition is not one quoted string')
        definitions.append(text.split('"')[1])
    index = {t: i for i, t in enumerate(vocabulary)}
    edges = sorted({(index[a], index[b]) for a, b, _ in graph.edges(keys=True) if a in index and b in index})
    return definitions, edges


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    common.add_argument(parser)
    root = common.data_root(parser.parse_args().data)
    from sentence_transformers import SentenceTransformer
    from torch_geometric.data import Data
    from torch_geometric.transforms import Compose, NormalizeFeatures, RandomLinkSplit, ToDevice
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    ontology = root / 'datasets' / 'ontology'
    vocabulary = json.loads((ontology / 'vocabulary.json').read_text())
    definitions, edges = definitions_and_edges(ontology / 'go-basic.obo', vocabulary)
    out = folder(root)
    out.mkdir(parents=True, exist_ok=True)
    encoder = SentenceTransformer(TEXT_MODEL, revision=TEXT_REVISION, device='cuda')
    text = torch.from_numpy(encoder.encode(definitions)).float()   # source: model.encode(definitions), defaults
    del encoder
    if text.shape != (len(vocabulary), LABEL_TEXT_DIM) or not torch.isfinite(text).all():
        raise ValueError(f'text features of shape {tuple(text.shape)} or not finite')
    edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()
    torch.save({'biobert': text, 'edge_index': edge_index}, out / 'graph.pt')

    seed_all(42)
    # `x` stands for the source's hierarchy matrix, which the biobert encoder never reads; a placeholder keeps the
    # transform's NormalizeFeatures step (which acts on x only) a no-op on the inputs that matter.
    data = Data(x=torch.zeros(len(vocabulary), 1), biobert=text, edge_index=edge_index, num_nodes=len(vocabulary))
    transform = Compose([NormalizeFeatures(), ToDevice('cuda'),
                         RandomLinkSplit(num_val=0.2, num_test=0.0, is_undirected=False, split_labels=True,
                                         add_negative_train_samples=True)])
    train, val, _ = transform(data)
    net = label_gae().cuda()
    optimizer = torch.optim.Adam(net.parameters(), lr=1e-4, weight_decay=1e-2)
    best, best_epoch, state = math.inf, None, None
    for epoch in range(500):
        net.train()
        optimizer.zero_grad()
        loss = net.recon_loss(net.encode(train.biobert, train.edge_index), train.pos_edge_label_index)
        loss.backward()
        optimizer.step()
        net.eval()
        with torch.no_grad():
            validation = float(net.recon_loss(net.encode(val.biobert, val.edge_index), val.pos_edge_label_index))
        if not (math.isfinite(validation) and math.isfinite(float(loss))):
            raise ValueError(f'non-finite loss at epoch {epoch}')
        if validation <= best:
            best, best_epoch = validation, epoch
            state = {k: t.detach().cpu().clone() for k, t in net.state_dict().items()}
        if epoch % 50 == 0:
            common.log(f'epoch {epoch}: train {float(loss):.4f}, validation {validation:.4f}')
    torch.save({'state_dict': state, 'epoch': best_epoch, 'valid_loss_min': best}, out / 'gcn.pt')
    common.log(f'label GCN: {len(vocabulary):,} terms, {len(edges):,} edges, best epoch {best_epoch} ({best:.4f})')


if __name__ == '__main__':
    main()
