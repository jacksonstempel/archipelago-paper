"""TransFew's full MF model and loss, copied from the official repository (MIT; LICENSE-TransFew).

Source: github.com/BioinfoMachineLearning/TransFew at commit 5a4c9d7d726c95891deb7f7789622749a44bf62a:
  - MLP, TFun_submodel, LabelEncoder, Attention, scatter_accross, TFun: models/model.py
  - the layer tables esm_layers_{freq,rare}_mf: models/config.py
  - get_weights (the loss weights): training.py

Changes, none of which alter the computation:
  - Shapes come from the arguments (number of outputs, frequent and rare columns) instead of module globals
    (official MF: 7,224 outputs; here the vocabulary and the frequent/rare split of the training set at hand).
  - The frozen label GCN and its graph are passed in instead of read from the authors' folders.
  - BP-only branches, unused sub-models (MSA, InterPro) and unused imports are dropped; torch_geometric is imported
    where the label GCN is built, so the rest of the module works without it.
  - TFun builds the frequent and rare sub-models in the official order, including its double construction, so a given
    torch seed gives the same initial weights as the official code.
"""
import math
from collections import OrderedDict

import torch
from torch import nn
from torch.nn import BatchNorm1d, GELU, LayerNorm, ReLU, Sigmoid, Softmax, Tanh

ESM_DIM = 5120
LABEL_TEXT_DIM = 768      # S-PubMedBert-MS-MARCO sentence embedding ("biobert" in the source)
LABEL_DIM = 1024          # LabelEncoder out_channels
ATTENTION_DIM = 256


def esm_layers(out):
    """models/config.py esm_layers_freq_mf / esm_layers_rare_mf with the last layer sized `out`:
    (in, out, bias, activation, norm, dropout, inputs)."""
    return [(ESM_DIM, 1024, True, 'gelu', 'batchnorm', 0.2, (0,)),
            (1024, 256, True, 'gelu', 'batchnorm', 0.2, (1,)),
            (256, out, True, 'gelu', 'batchnorm', 'none', (2,))]


def getactfn(actfn='relu'):
    return {'relu': ReLU, 'gelu': GELU, 'tanh': Tanh}[actfn]()


def getnorm(norm='layernorm', norm_shape=0):
    return {'layernorm': LayerNorm, 'batchnorm': BatchNorm1d}[norm](norm_shape)


class MLP(nn.Module):
    """models/model.py MLP (aggregate=None path)."""

    def __init__(self, layers, in_dropout=0.0):
        super().__init__()
        output_dim, input_dim = layers[-1][1], layers[0][0]
        self.in_dropout = nn.Dropout(in_dropout)
        self.layers = nn.ModuleList()
        self.layer_inputs = []
        for pos, d in enumerate(layers):
            tmp = OrderedDict()
            tmp[f'mlp_{pos}'] = nn.Linear(d[0], d[1], bias=d[2])
            if d[4] != 'none':
                tmp[f'norm_{pos}'] = getnorm(d[4], d[1])
            if d[3] != 'none':
                tmp[f'act_fn_{pos}'] = getactfn(d[3])
            if d[5] != 'none':
                tmp[f'dropout_{pos}'] = nn.Dropout(p=d[5])
            self.layer_inputs.append(d[6])
            self.layers.append(nn.Sequential(tmp))
        self.skip_1 = nn.Linear(input_dim, output_dim)
        self.skip_final = nn.Linear(output_dim * 2, output_dim)

    def forward(self, features, residual=True):
        features = self.in_dropout(features)
        features_array = [features]
        for lay, inputs in zip(self.layers, self.layer_inputs):
            x = lay(torch.cat([features_array[i] for i in inputs], 1))
            features_array.append(x)
        if residual:
            x = self.skip_final(torch.cat([x, self.skip_1(features_array[0])], 1))
        return x


class TFun_submodel(nn.Module):
    """models/model.py TFun_submodel, sub_model 'esm2_t48', full=True (no sigmoid)."""

    def __init__(self, out):
        super().__init__()
        self.esm_mlp = MLP(esm_layers(out), in_dropout=0.1)

    def forward(self, data):
        return self.esm_mlp(data, residual=True)


class LabelEncoder(nn.Module):
    """models/model.py LabelEncoder (features='biobert')."""

    def __init__(self, in_channels, out_channels):
        from torch_geometric.nn import GCNConv
        super().__init__()
        self.conv1 = GCNConv(in_channels, 2 * out_channels, cached=True)
        self.conv2 = GCNConv(2 * out_channels, out_channels, cached=True)
        self.dropout = nn.Dropout(0.1)
        self.relu = nn.ReLU()
        self.tanh = nn.Tanh()
        self.bn = BatchNorm1d(2 * out_channels)

    def forward(self, x, edge_index):
        x = self.dropout(x)
        x = self.conv1(x, edge_index)
        x = self.bn(x)
        x = self.tanh(x)
        return self.conv2(x, edge_index)


def label_gae():
    from torch_geometric.nn import GAE
    return GAE(LabelEncoder(LABEL_TEXT_DIM, LABEL_DIM))


class Attention(nn.Module):
    """models/model.py Attention: proteins as queries, GO terms as keys and values."""

    def __init__(self, concat_hidden, label_dimension, embedding_size, out_proj):
        super().__init__()
        self.q = nn.Linear(concat_hidden, embedding_size)
        self.k = nn.Linear(label_dimension, embedding_size)
        self.v = nn.Linear(label_dimension, embedding_size)
        self.final = nn.Linear(embedding_size, out_proj)
        self.softmax = Softmax(dim=-1)

    def forward(self, x_1, x_2):
        q, k, v = self.q(x_1), self.k(x_2), self.v(x_2)
        scores = torch.matmul(q, k.transpose(-2, -1)) / math.sqrt(k.size()[-1])
        values = torch.matmul(self.softmax(scores), v)
        return self.final(values + q)


def scatter_accross(shape, freq_scores, rare_scores, freq_indicies, rare_indicies):
    """models/model.py scatter_accross (MF: no rare_2 group)."""
    batch = freq_scores.shape[0]
    des = torch.zeros(batch, shape, dtype=freq_scores.dtype, device=freq_scores.device)
    des = des.scatter_(1, freq_indicies.repeat(batch, 1), freq_scores)
    return des.scatter_(1, rare_indicies.repeat(batch, 1), rare_scores)


class TFun(nn.Module):
    """models/model.py TFun for MF. `freq` and `rare` are the output columns of the two branches (together every
    column once); `label_state` is the pretrained GAE state dict; `graph` holds `biobert` (terms x 768, vocabulary
    order) and `edge_index`. The caller freezes the label GAE's parameters, as training.py does."""

    def __init__(self, n_terms, freq, rare, label_state, graph):
        super().__init__()
        if sorted(torch.cat([freq, rare]).tolist()) != list(range(n_terms)):
            raise ValueError('frequent and rare columns must cover every output exactly once')
        self.out_shape = n_terms
        self.register_buffer('freq_indicies', freq.clone())
        self.register_buffer('rare_indicies', rare.clone())
        self.dropout = nn.Dropout(0.1)    # unused in forward, as in the source
        self.dropout2 = nn.Dropout(0.1)
        self.esm_freq_mlp = TFun_submodel(len(freq))   # official: load_submodel(...) then TFun_submodel(...)
        self.esm_freq_mlp = TFun_submodel(len(freq))
        self.esm_rare_mlp = TFun_submodel(len(rare))   # official: TFun_submodel(...) then load_submodel(...)
        self.esm_rare_mlp = TFun_submodel(len(rare))
        self.label_embedding = label_gae()
        self.label_embedding.load_state_dict(label_state)
        self.joint_embedding = Attention(concat_hidden=n_terms, label_dimension=LABEL_DIM,
                                         embedding_size=ATTENTION_DIM, out_proj=n_terms)
        self.register_buffer('label_text', graph['biobert'].float().clone(), persistent=False)
        self.register_buffer('label_edges', graph['edge_index'].clone(), persistent=False)
        self.sigmoid = Sigmoid()

    def forward(self, x):
        scores = scatter_accross(self.out_shape, self.esm_freq_mlp(x), self.esm_rare_mlp(x),
                                 self.freq_indicies, self.rare_indicies)
        z = self.label_embedding.encode(self.label_text, self.label_edges)
        return self.sigmoid(self.joint_embedding(scores, z))


def get_weights(labels, alpha):
    """training.py get_weights: log2(sum of batch label counts / (count + 0.01)) ** alpha, per column."""
    label_frequencies = torch.sum(labels, dim=0)
    total_count = torch.sum(label_frequencies, dim=0)
    label_frequencies = torch.add(label_frequencies, 1e-2)
    return torch.pow(torch.log2(torch.div(total_count, label_frequencies)), alpha)


def weighted_bce(output, labels, alpha=3):
    """training.py loss of the full model: (BCELoss(reduction='none')(output, labels) * get_weights(labels, 3)).mean()."""
    return (torch.nn.functional.binary_cross_entropy(output, labels, reduction='none') * get_weights(labels, alpha)).mean()
