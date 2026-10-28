"""DeepGO-SE's model and axiom loader, copied unchanged from the official code (BSD 3-Clause, LICENSE-deepgo2 here).

Source: github.com/bio-ontology-research-group/deepgo2, commit 6c77c985ef9453548dca629830c044b2dc1b4982:
- Residual, MLPBlock, BaseDeepGOModel: deepgo/base.py (whole file).
- DeepGOModel: deepgo/models.py (only this class; the file's GAT model needs DGL and is not used).
- load_normal_forms: deepgo/data.py (only this function).
Only docstrings were shortened; code, defaults and quirks are the official ones (for example nf4_loss reads go_rad with the
relation index). The training step in engine.py computes exactly this model's loss and gradients.
"""
import math

import torch as th
from torch import nn


class Residual(nn.Module):
    """A residual layer that adds the output of a function to its input."""

    def __init__(self, fn):
        super().__init__()
        self.fn = fn

    def forward(self, x):
        return x + self.fn(x)


class MLPBlock(nn.Module):
    """One fully connected layer with activation, layer norm and dropout."""
    def __init__(self, in_features, out_features, bias=True, layer_norm=True, dropout=0.1, activation=nn.ReLU):
        super().__init__()
        self.linear = nn.Linear(in_features, out_features, bias)
        self.activation = activation()
        self.layer_norm = nn.LayerNorm(out_features) if layer_norm else None
        self.dropout = nn.Dropout(dropout) if dropout else None

    def forward(self, x):
        x = self.activation(self.linear(x))
        if self.layer_norm:
            x = self.layer_norm(x)
        if self.dropout:
            x = self.dropout(x)
        return x


class BaseDeepGOModel(nn.Module):
    """A base DeepGO model with ElEmbeddings loss functions."""

    def __init__(self, input_length, nb_gos, nb_zero_gos, nb_rels, device, hidden_dim=2560, embed_dim=2560, margin=0.1):
        super().__init__()
        self.nb_gos = nb_gos
        self.nb_zero_gos = nb_zero_gos
        self.nb_rels = nb_rels
        # ELEmbedding Model Layers
        self.embed_dim = embed_dim
        # Create additional index for hasFunction relation
        self.hasFuncIndex = th.LongTensor([nb_rels]).to(device)
        # Embedding layer for all classes in GO
        self.go_embed = nn.Embedding(nb_gos + nb_zero_gos, embed_dim)
        self.go_norm = nn.BatchNorm1d(embed_dim)
        # Initialize embedding layers
        k = math.sqrt(1 / embed_dim)
        nn.init.uniform_(self.go_embed.weight, -k, k)
        self.go_rad = nn.Embedding(nb_gos + nb_zero_gos, 1)
        nn.init.uniform_(self.go_rad.weight, -k, k)
        self.rel_embed = nn.Embedding(nb_rels + 1, embed_dim)
        nn.init.uniform_(self.rel_embed.weight, -k, k)
        # indices of all GO classes
        self.all_gos = th.arange(self.nb_gos).to(device)
        self.margin = margin

    def forward(self, features):
        raise NotImplementedError

    def el_loss(self, go_normal_forms):
        nf1, nf2, nf3, nf4 = go_normal_forms
        loss = self.nf1_loss(nf1)
        if len(nf2):
            loss += self.nf2_loss(nf2)
        if len(nf3):
            loss += self.nf3_loss(nf3)
        if len(nf4):
            loss += self.nf4_loss(nf4)
        return loss

    def class_dist(self, data):
        c = self.go_norm(self.go_embed(data[:, 0]))
        d = self.go_norm(self.go_embed(data[:, 1]))
        rc = th.abs(self.go_rad(data[:, 0]))
        rd = th.abs(self.go_rad(data[:, 1]))
        dist = th.linalg.norm(c - d, dim=1, keepdim=True) + rc - rd
        return dist

    def nf1_loss(self, data):
        """C subclassOf D"""
        pos_dist = self.class_dist(data)
        loss = th.mean(th.relu(pos_dist - self.margin))
        return loss

    def nf2_loss(self, data):
        """C and D subclassOf E"""
        c = self.go_norm(self.go_embed(data[:, 0]))
        d = self.go_norm(self.go_embed(data[:, 1]))
        e = self.go_norm(self.go_embed(data[:, 2]))
        rc = th.abs(self.go_rad(data[:, 0]))
        rd = th.abs(self.go_rad(data[:, 1]))
        re = th.abs(self.go_rad(data[:, 2]))

        sr = rc + rd
        dst = th.linalg.norm(c - d, dim=1, keepdim=True)
        dst2 = th.linalg.norm(e - c, dim=1, keepdim=True)
        dst3 = th.linalg.norm(e - d, dim=1, keepdim=True)
        loss = th.mean(th.relu(dst - sr - self.margin)
                       + th.relu(dst2 - rc - self.margin)
                       + th.relu(dst3 - rd - self.margin))
        return loss

    def nf3_loss(self, data):
        """R some C subClassOf D"""
        rE = self.rel_embed(data[:, 0])
        c = self.go_norm(self.go_embed(data[:, 1]))
        d = self.go_norm(self.go_embed(data[:, 2]))
        rc = th.abs(self.go_rad(data[:, 1]))
        rd = th.abs(self.go_rad(data[:, 2]))

        rSomeC = c + rE
        euc = th.linalg.norm(rSomeC - d, dim=1, keepdim=True)
        loss = th.mean(th.relu(euc + rc - rd - self.margin))
        return loss

    def nf4_loss(self, data):
        """C subclassOf R some D"""
        c = self.go_norm(self.go_embed(data[:, 0]))
        rE = self.rel_embed(data[:, 1])
        d = self.go_norm(self.go_embed(data[:, 2]))

        rc = th.abs(self.go_rad(data[:, 1]))
        rd = th.abs(self.go_rad(data[:, 2]))
        sr = rc + rd
        # c should intersect with d + r
        rSomeD = d + rE
        dst = th.linalg.norm(c - rSomeD, dim=1, keepdim=True)
        loss = th.mean(th.relu(dst - sr - self.margin))
        return loss


class DeepGOModel(BaseDeepGOModel):
    """DeepGO model with ElEmbeddings loss functions."""

    def __init__(self, input_length, nb_gos, nb_zero_gos, nb_rels, device, hidden_dim=2560, embed_dim=2560, margin=0.1):
        super().__init__(input_length, nb_gos, nb_zero_gos, nb_rels, device, hidden_dim, embed_dim, margin)
        # MLP Layers to project the input protein
        net = []
        net.append(MLPBlock(input_length, embed_dim))
        net.append(Residual(MLPBlock(embed_dim, embed_dim)))
        self.net = nn.Sequential(*net)

    def forward(self, features):
        x = self.net(features)
        go_embed = self.go_embed(self.all_gos)
        hasFunc = self.rel_embed(self.hasFuncIndex)
        hasFuncGO = go_embed + hasFunc
        go_rad = th.abs(self.go_rad(self.all_gos).view(1, -1))
        x = th.matmul(x, hasFuncGO.T) + go_rad
        logits = th.sigmoid(x)
        return logits


def load_normal_forms(go_file, terms_dict):
    """Parses the normalized axioms (axioms/Normalizer.groovy output). Classes of terms_dict get its indices; every other
    class gets the next free index (zero classes). Lines that fit no branch (e.g. NF1 with a non-10-character class) are
    skipped, as in the official loader."""
    nf1 = []
    nf2 = []
    nf3 = []
    nf4 = []
    relations = {}
    zclasses = {}

    def get_index(go_id):
        if go_id in terms_dict:
            index = terms_dict[go_id]
        elif go_id in zclasses:
            index = zclasses[go_id]
        else:
            zclasses[go_id] = len(terms_dict) + len(zclasses)
            index = zclasses[go_id]
        return index

    def get_rel_index(rel_id):
        if rel_id not in relations:
            relations[rel_id] = len(relations)
        return relations[rel_id]

    with open(go_file) as f:
        for line in f:
            line = line.strip().replace('_', ':')
            if line.find('SubClassOf') == -1:
                continue
            left, right = line.split(' SubClassOf ')
            # C SubClassOf D
            if len(left) == 10 and len(right) == 10:
                go1, go2 = left, right
                nf1.append((get_index(go1), get_index(go2)))
            elif left.find('and') != -1:  # C and D SubClassOf E
                go1, go2 = left.split(' and ')
                go3 = right
                nf2.append((get_index(go1), get_index(go2), get_index(go3)))
            elif left.find('some') != -1:  # R some C SubClassOf D
                rel, go1 = left.split(' some ')
                go2 = right
                nf3.append((get_rel_index(rel), get_index(go1), get_index(go2)))
            elif right.find('some') != -1:  # C SubClassOf R some D
                go1 = left
                rel, go2 = right.split(' some ')
                nf4.append((get_index(go1), get_rel_index(rel), get_index(go2)))
    return nf1, nf2, nf3, nf4, relations, zclasses
