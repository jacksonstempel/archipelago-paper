"""The released DPFunc model and loss, copied without change of computation.

Source: github.com/CSUBioGroup/DPFunc at commit c5faf4df4e29b7d264ccb14d88254a136b54f4b2, MIT License, Copyright (c) 2024
CSUBioGroup (LICENSE-DPFunc in this folder): `DPFunc/models.py` (inter_model, transformer_block, GCN,
combine_inter_model) and `DPFunc/model_utils.py` (FocalLoss). Unused imports, comments and blank lines removed; nothing
else changed. The one runtime adaptation (residues without contacts) is applied by train.py.
"""
import dgl
import dgl.nn as dglnn
import torch
import torch.nn as nn
import torch.nn.functional as F


class inter_model(nn.Module):
    def __init__(self, input_size, hidden_size):
        super(inter_model, self).__init__()
        self.embedding_layer = nn.EmbeddingBag(input_size, hidden_size, mode='sum', include_last_offset=True)
        self.linearLayer = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.Dropout(0.3),
            nn.ReLU(),
            nn.Linear(hidden_size, hidden_size),
            nn.Dropout(0.3),
            nn.ReLU()
        )

    def forward(self, inter_feature):
        inter_feature = F.relu(self.embedding_layer(*inter_feature))
        inter_feature = self.linearLayer(inter_feature)
        return inter_feature


class transformer_block(nn.Module):
    def __init__(self, in_dim, hidden_dim, head=1):
        super(transformer_block, self).__init__()
        self.head = head
        self.trans_q_list = nn.ModuleList([nn.Linear(in_dim, hidden_dim, bias=False) for _ in range(head)])
        self.trans_k_list = nn.ModuleList([nn.Linear(in_dim, hidden_dim, bias=False) for _ in range(head)])
        self.trans_v_list = nn.ModuleList([nn.Linear(in_dim, hidden_dim, bias=False) for _ in range(head)])
        self.concat_trans = nn.Linear((hidden_dim)*head, hidden_dim, bias=False)
        self.ff = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim*2),
            nn.ReLU(),
            nn.Linear(hidden_dim*2, hidden_dim)
        )
        self.layernorm = nn.LayerNorm(in_dim)

    def forward(self, g, residue_h, inter_h):
        multi_output = []
        for i in range(self.head):
            q = self.trans_q_list[i](residue_h)
            k = self.trans_k_list[i](inter_h)
            v = self.trans_v_list[i](residue_h)
            att = torch.sum(torch.mul(q, k)/torch.sqrt(torch.tensor(1280.0)), dim=1, keepdim=True)
            with g.local_scope():
                g.ndata['att'] = att.reshape(-1)
                alpha = dgl.softmax_nodes(g, 'att').reshape((v.size(0), 1))
                tp = v * alpha
            multi_output.append(tp)
        multi_output = torch.cat(multi_output, dim=1)
        multi_output = self.concat_trans(multi_output)
        multi_output = self.layernorm(multi_output + residue_h)
        multi_output = self.layernorm(self.ff(multi_output)+multi_output)
        return multi_output


class GCN(nn.Module):
    def __init__(self, in_dim, hidden_dim, n_classes, head):
        super(GCN, self).__init__()
        self.dropout = nn.Dropout(0.3)
        self.bn1 = nn.BatchNorm1d(hidden_dim)
        self.bn2 = nn.BatchNorm1d(hidden_dim)
        self.conv1 = dglnn.GraphConv(in_dim, hidden_dim)
        self.conv2 = dglnn.GraphConv(hidden_dim, hidden_dim)
        self.transformer_block = transformer_block(hidden_dim, hidden_dim, head)

    def forward(self, g, h, inter_f):
        with g.local_scope():
            g.ndata['h'] = h
            init_avg_h = dgl.mean_nodes(g, 'h')
        pre = h
        h = self.bn1(h)
        h = pre + self.dropout(F.relu(self.conv1(g, h)))
        pre = h
        h = self.bn2(h)
        h = pre + self.dropout(F.relu(self.conv2(g, h)))
        with g.local_scope():
            g.ndata['inter'] = dgl.broadcast_nodes(g, inter_f)
            residue_h = h
            inter_h = g.ndata['inter']
            hg = self.transformer_block(g, residue_h, inter_h)
            g.ndata['output'] = hg
            readout = dgl.sum_nodes(g, "output")
            return readout, init_avg_h


class combine_inter_model(nn.Module):
    def __init__(self, inter_size, inter_hid, graph_size, graph_hid, label_num, head):
        super(combine_inter_model, self).__init__()
        self.inter_embedding = inter_model(inter_size, inter_hid)
        self.GNN = GCN(graph_size, graph_hid, label_num, head)
        self.classify = nn.Sequential(
            nn.BatchNorm1d(graph_size+graph_hid),
            nn.Linear(graph_size+graph_hid, (graph_size+graph_hid)*2),
            nn.Dropout(0.3),
            nn.ReLU(),
            nn.Linear((graph_size+graph_hid)*2, (graph_size+graph_hid)*2),
            nn.Dropout(0.3),
            nn.ReLU(),
            nn.Linear((graph_size+graph_hid)*2, label_num)
        )

    def forward(self, inter_feature, graph, graph_h):
        inter_feature = self.inter_embedding(inter_feature)
        graph_feature, init_feature = self.GNN(graph, graph_h, inter_feature)
        return self.classify(torch.cat((init_feature, graph_feature), 1))


class FocalLoss(torch.nn.Module):
    def __init__(self, gamma=2, alpha=None):
        super(FocalLoss, self).__init__()
        self.gamma = gamma
        self.alpha = alpha

    def forward(self, pred, target):
        ce_loss = F.binary_cross_entropy_with_logits(pred, target, reduction='none')
        pt = torch.exp(-ce_loss)
        focal_loss = ((1 - pt) ** self.gamma * ce_loss).mean()
        if self.alpha is not None:
            alpha_t = self.alpha[target]
            focal_loss = alpha_t * focal_loss
        return focal_loss
