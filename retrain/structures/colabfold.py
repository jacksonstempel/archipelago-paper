"""Structure bank, part 2: AlphaFold2 (ColabFold 1.6.3) for every protein without an AlphaFold DB model of its sequence.

    python -m retrain.structures.colabfold msa --databases DIR [--mmseqs BIN] [--threads 32] [--part i --parts n]
    python -m retrain.structures.colabfold fold [--part i --parts n]

The targets are needs_colabfold.tsv (written by retrain.structures.afdb fetch): one target per distinct modeled sequence
(residues 1-1,400 above 2,700 residues; U/O/B/Z/X replaced by C/K/D/E/A), named by its SHA-256.

msa   ColabFold 1.6.3's own MMseqs2 search (colabfold.mmseqs.search, as colabfold_search runs it) against UniRef30 2302
      and colabfold_envdb_202108 (the databases of ColabFold's setup_databases.sh): environmental database on, no
      templates, filter on, prefilter mode 0, unpaired -> colabfold/msas/<target>.a3m. As in the paper, ColabFold's
      quoted k-score ('seq:96,prof:80', which MMseqs2 18 would otherwise ignore) is passed unquoted, and MMseqs2 must be
      release 18 (8cc5c). CPU.
fold  colabfold_batch on the MSAs: alphafold2_ptm, model 1 only, one seed (0), 3 recycles, no templates, no relaxation,
      recompile padding 128 -> colabfold/models/<target[:2]>/<target>.pdb.gz. One GPU per part.

Both stages skip finished targets, and --part i --parts n splits the targets (sorted by length) across jobs; a fold job
that stopped keeps its work folder, and colabfold_batch skips the folds it finished. Needs ColabFold 1.6.3
(`pip install "colabfold[alphafold]==1.6.3"` with a CUDA JAX) and MMseqs2 release 18. ColabFold is MIT licensed; the
AlphaFold2 parameters (downloaded by ColabFold on first use) are CC BY 4.0.

The paper's bank used this recipe but, to save compute, reused folds made earlier in the project from the same MSAs:
folds with the standard kernels, 1,762 with ColabFold's fused (Pallas) kernels and 34,095 with the fused kernels and
bucket padding (a fast path whose folds agree with standard-kernel folds at 97.3% Foldseek 3Di identity); the last ~944
targets were also folded with that fast path. A bank rebuilt here uses the standard kernels for every target, so some
structures differ slightly from the paper's.
"""
import argparse
import gzip
import hashlib
import shutil
import subprocess
import sys
import time

from retrain.common import add_argument, data_root, log
from retrain.structures.afdb import load_requests
from retrain.structures.common import bank, read_structure, read_tsv, target_sequence, write_bytes

SUFFIX = '_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb'
FOLD = ['--model-type', 'alphafold2_ptm', '--num-models', '1', '--model-order', '1', '--num-seeds', '1',
        '--random-seed', '0', '--num-recycle', '3', '--num-relax', '0', '--recompile-padding', '128']
MMSEQS_RELEASE = '8cc5ce367b5638c4306c2d7cfc652dd099a4643f'  # release 18


def targets(root, part, parts):
    """{target: modeled sequence} of this part (targets sorted by length, then name; every n-th from i)."""
    sequences = {r['seq_id']: r['sequence'] for r in load_requests(root)}
    found = {}
    for r in read_tsv(bank(root) / 'needs_colabfold.tsv'):
        target = target_sequence(sequences[r['seq_id']])
        if hashlib.sha256(target.encode()).hexdigest() != r['target']:
            raise ValueError(f"needs_colabfold.tsv: {r['seq_id']} names another target; run the afdb fetch stage again")
        found[r['target']] = target
    ordered = sorted(found, key=lambda t: (len(found[t]), t))
    return {t: found[t] for t in ordered[part::parts]}


def msa_path(root, target):
    return bank(root) / 'colabfold' / 'msas' / f'{target}.a3m'


def model_path(root, target):
    return bank(root) / 'colabfold' / 'models' / target[:2] / f'{target}.pdb.gz'


def query_sequence(a3m):
    """The first sequence of an A3M file (ColabFold's query)."""
    rows, seen = [], False
    for row in open(a3m):
        if row.startswith('>'):
            if seen:
                break
            seen = True
        elif seen and not row.startswith('#'):
            rows.append(row.strip())
    return ''.join(rows)


def msa(root, arguments):
    todo = {t: s for t, s in targets(root, arguments.part, arguments.parts).items() if not msa_path(root, t).exists()}
    if not todo:
        return log('msa: nothing to do')
    work = bank(root) / 'colabfold' / f'search_{arguments.part}'
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    (work / 'queries.fasta').write_text(''.join(f'>{t}\n{s}\n' for t, s in todo.items()))
    version = subprocess.check_output([arguments.mmseqs, 'version'], text=True).strip()
    if version != MMSEQS_RELEASE:
        raise SystemExit(f'MMseqs2 {version} is not release 18 ({MMSEQS_RELEASE}), which made the paper\'s MSAs')
    log(f'msa: {len(todo):,} targets, {sum(map(len, todo.values())):,} residues')
    started = time.monotonic()
    subprocess.run([sys.executable, '-m', 'retrain.structures.colabfold', 'search-child', work / 'queries.fasta',
                    arguments.databases, work / 'out', '--mmseqs', arguments.mmseqs, '--threads', str(arguments.threads),
                    '--use-env', '1', '--use-templates', '0', '--filter', '1', '--prefilter-mode', '0',
                    '--db-load-mode', '0', '--pair-mode', 'unpaired'], check=True)
    for t, s in todo.items():
        a3m = work / 'out' / f'{t}.a3m'
        if query_sequence(a3m) != s:
            raise ValueError(f'{a3m}: the MSA query is not the target sequence')
        msa_path(root, t).parent.mkdir(parents=True, exist_ok=True)
        a3m.replace(msa_path(root, t))
    shutil.rmtree(work)
    log(f'msa: done in {(time.monotonic() - started) / 3600:.2f} h')


def fold(root, arguments):
    todo = {t: s for t, s in targets(root, arguments.part, arguments.parts).items() if not model_path(root, t).exists()}
    missing = [t for t in todo if not msa_path(root, t).exists()]
    if missing:
        raise SystemExit(f'{len(missing):,} targets have no MSA yet; run the msa stage first')
    if not todo:
        return log('fold: nothing to do')
    work = bank(root) / 'colabfold' / f'fold_{arguments.part}'  # kept after a crash: colabfold_batch skips finished folds
    (work / 'msas').mkdir(parents=True, exist_ok=True)
    for t in todo:
        if not (work / 'msas' / f'{t}.a3m').exists():
            shutil.copyfile(msa_path(root, t), work / 'msas' / f'{t}.a3m')
    log(f'fold: {len(todo):,} targets, {sum(map(len, todo.values())):,} residues')
    started = time.monotonic()
    subprocess.run(['colabfold_batch', *FOLD, work / 'msas', work / 'out'], check=True)
    for t, s in todo.items():
        payload = (work / 'out' / f'{t}{SUFFIX}').read_bytes()
        sequence, _, _ = read_structure(payload)
        if sequence != s:
            raise ValueError(f'the fold of {t} has another sequence')
        write_bytes(model_path(root, t), gzip.compress(payload, 6))
    shutil.rmtree(work)
    log(f'fold: {len(todo):,} targets in {(time.monotonic() - started) / 3600:.2f} h')


def search_child(argv):
    """ColabFold 1.6.3's own search (colabfold.mmseqs.search), as the paper ran it: ColabFold passes the k-score with
    literal shell quotes ('seq:96,prof:80'), which MMseqs2 18 ignores, so the quotes are removed and the k-score applies;
    searches get --split-memory-limit 160G."""
    import colabfold.mmseqs.search as native
    original = native.run_mmseqs

    def bounded(binary, params):
        params = list(params)
        if '--k-score' in params:
            i = params.index('--k-score') + 1
            if params[i] != "'seq:96,prof:80'":
                raise RuntimeError(f'unexpected k-score {params[i]}: another ColabFold version?')
            params[i] = 'seq:96,prof:80'
        if params[0] == 'search':
            params += ['--split-memory-limit', '160G']
        return original(binary, params)

    native.run_mmseqs = bounded
    sys.argv = ['colabfold_search', *argv]
    native.main()


def main():
    if sys.argv[1:2] == ['search-child']:  # the msa stage's subprocess
        return search_child(sys.argv[2:])
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('stage', choices=['msa', 'fold'])
    parser.add_argument('--databases', help='msa: folder of the ColabFold MMseqs2 databases (setup_databases.sh)')
    parser.add_argument('--threads', type=int, default=32, help='msa: MMseqs2 threads')
    parser.add_argument('--mmseqs', default='mmseqs', help='msa: the MMseqs2 release 18 binary')
    parser.add_argument('--part', type=int, default=0)
    parser.add_argument('--parts', type=int, default=1)
    add_argument(parser)
    arguments = parser.parse_args()
    if arguments.stage == 'msa' and not arguments.databases:
        parser.error('msa needs --databases')
    globals()[arguments.stage](data_root(arguments.data), arguments)


if __name__ == '__main__':
    main()
