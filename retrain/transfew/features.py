"""TransFew's input features: ESM2-15B layer 48, mean over residues of non-overlapping 1,021-residue chunks, FP32.

The extractor of the official DataGen/Embeddings.py and predict.py (crop each sequence into consecutive chunks of at
most 1,021 residues, run each chunk alone, concatenate the residue outputs, mean of layer 48): Transformers EsmModel
with the official weights facebook/esm2_t48_15B_UR50D (revision 5fbca39631164edc1d402a5aa369f982f72ee282), FP32
weights and inference with TF32 off, `last_hidden_state` (layer 48 after the final layer norm, as fair-esm's
representation 48). Sequences are passed as they are (the tokenizer maps unknown residues).

Every distinct sequence of the named datasets is embedded once. Stages:
  embed   one GPU worker: shard K of N of the sequences (run N workers in parallel, e.g. one per GPU); work is saved in
          chunks of 500 proteins under inputs/transfew/esm2_15b_l48/parts/, so a stopped worker resumes
  merge   one matrix of every sequence: inputs/transfew/esm2_15b_l48/vectors.npy + vectors.ids.txt (archipelago.io)

    python -m retrain.transfew.features embed --shard 0 --shards 4 [--data DIR] [--sets ...]
    python -m retrain.transfew.features merge [--data DIR] [--sets ...]

Speed: about 0.3 to 1 s per protein on one H100 (FP32, one protein at a time; the model needs about 60 GB), so the
default sets (about 310,000 distinct sequences) take on the order of 40 GPU-hours.
"""
import argparse
import time
from pathlib import Path

import numpy as np

from archipelago.io import write_matrix
from retrain import common

MODEL = 'facebook/esm2_t48_15B_UR50D'
REVISION = '5fbca39631164edc1d402a5aa369f982f72ee282'
CHUNK = 1021
DIMENSION = 5120
SETS = ('original_transfew', 'goa', 'hq', 'archipelago', 'development', 'test')
PART = 500  # proteins per saved chunk


def folder(root):
    return Path(root) / 'inputs' / 'transfew' / 'esm2_15b_l48'


def sequences(root, sets):
    """{seq_id: sequence} of every protein of the datasets, sorted by seq_id."""
    found = {}
    for name in sets:
        ds = common.load(root, name, labels=False)
        found.update(zip(ds.seq_ids, ds.sequences))
    return dict(sorted(found.items()))


class Extractor:
    def __init__(self):
        import torch
        from transformers import AutoTokenizer, EsmModel
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        self.torch = torch
        self.tokenizer = AutoTokenizer.from_pretrained(MODEL, revision=REVISION)
        self.model = EsmModel.from_pretrained(MODEL, revision=REVISION, torch_dtype=torch.float32,
                                              low_cpu_mem_usage=True, add_pooling_layer=False).eval().cuda()

    def __call__(self, sequence):
        torch = self.torch
        total = torch.zeros(DIMENSION, device='cuda', dtype=torch.float32)
        with torch.no_grad():
            for start in range(0, len(sequence), CHUNK):
                chunk = sequence[start:start + CHUNK]
                tokens = self.tokenizer(chunk, return_tensors='pt')
                if tokens['input_ids'].shape != (1, len(chunk) + 2):
                    raise ValueError(f'unexpected tokenization of a {len(chunk)}-residue chunk')
                hidden = self.model(**{k: v.cuda() for k, v in tokens.items()}).last_hidden_state
                total += hidden[0, 1:len(chunk) + 1].sum(0)
        vector = (total / len(sequence)).cpu().numpy()
        if not np.isfinite(vector).all():
            raise ValueError('non-finite embedding')
        return vector


def embed(root, sets, shard, shards):
    todo = list(sequences(root, sets).items())[shard::shards]
    parts = folder(root) / 'parts' / f'{shard}-of-{shards}'
    parts.mkdir(parents=True, exist_ok=True)
    pending = [a for a in range(0, len(todo), PART) if not (parts / f'{a:08d}.npy').exists()]
    common.log(f'shard {shard} of {shards}: {len(todo):,} proteins, {len(pending)} of '
               f'{-(-len(todo) // PART)} chunks to do')
    if not pending:
        return
    extract = Extractor()
    started, residues = time.monotonic(), 0
    left = sum(len(s) for a in pending for _, s in todo[a:a + PART])
    for a in pending:
        block = todo[a:a + PART]
        vectors = np.stack([extract(s) for _, s in block])
        (parts / f'{a:08d}.ids.txt').write_text(''.join(f'{h}\n' for h, _ in block))
        np.save(parts / f'{a:08d}.tmp.npy', vectors)
        (parts / f'{a:08d}.tmp.npy').replace(parts / f'{a:08d}.npy')  # last: a chunk counts as done only now
        residues += sum(len(s) for _, s in block)
        rate = residues / (time.monotonic() - started)
        common.log(f'chunk {a // PART}: {rate:,.0f} residues/s, about {(left - residues) / rate / 3600:.1f} h left')


def merge(root, sets):
    wanted = sequences(root, sets)
    ids, blocks = [], []
    for shard_folder in sorted((folder(root) / 'parts').iterdir()):
        for path in sorted(shard_folder.glob('*[0-9].npy')):
            block = np.load(path)
            names = path.with_suffix('.ids.txt').read_text().split()
            if block.shape != (len(names), DIMENSION):
                raise ValueError(f'{path}: shape {block.shape} for {len(names)} identifiers')
            ids += names
            blocks.append(block)
    extra = set(ids) - set(wanted)
    if extra:
        raise ValueError(f'{len(extra):,} embedded sequences belong to none of the datasets; were the parts made for '
                         'other --sets?')
    missing = set(wanted) - set(ids)
    if missing or len(set(ids)) != len(ids):
        raise ValueError(f'{len(missing):,} sequences have no vector, or a sequence was embedded twice; run embed for '
                         'every shard with the same --sets and --shards')
    write_matrix(folder(root) / 'vectors', ids, np.concatenate(blocks).astype(np.float32))
    common.log(f'wrote {folder(root) / "vectors.npy"}: {len(ids):,} proteins')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('stage', choices=['embed', 'merge'])
    common.add_argument(parser)
    parser.add_argument('--sets', nargs='+', default=list(SETS), help='datasets whose proteins are embedded')
    parser.add_argument('--shard', type=int, default=0)
    parser.add_argument('--shards', type=int, default=1)
    args = parser.parse_args()
    root = common.data_root(args.data)
    if args.stage == 'embed':
        embed(root, args.sets, args.shard, args.shards)
    else:
        merge(root, args.sets)


if __name__ == '__main__':
    main()
