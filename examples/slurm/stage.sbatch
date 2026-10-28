#!/bin/bash
# Example only: run one pipeline stage as a Slurm job. Every stage also runs without Slurm, as the same command.
# Add your cluster's partition, account and QoS lines, and size the resources to the stage (README.md lists them).
#
#   sbatch --export=ALL,ARCHIPELAGO_DATA=/path/to/data examples/slurm/stage.sbatch retrain.finch.train --training-set archipelago --seed 9301
#   sbatch --export=ALL,ARCHIPELAGO_DATA=/path/to/data examples/slurm/stage.sbatch scripts/embed.py --set hq --encoder esm2_15b_l48
#
# Several fits as one array job, e.g. DeepGO-SE's ten worlds of one training set (task k trains world k):
#
#   sbatch --array=0-9 --export=ALL,ARCHIPELAGO_DATA=/path/to/data examples/slurm/stage.sbatch \
#       retrain.deepgose.train world archipelago '$SLURM_ARRAY_TASK_ID'
#
#SBATCH --job-name=archipelago
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --output=archipelago-%x-%A_%a.log
set -euo pipefail
: "${ARCHIPELAGO_DATA:?set ARCHIPELAGO_DATA to the data root}"
args=()
for a in "${@:2}"; do args+=("$(eval echo "$a")"); done  # expands a quoted '$SLURM_ARRAY_TASK_ID'
echo "start $(date) on $(hostname): $1 ${args[*]}"
if [[ $1 == *.py ]]; then python -u "$1" "${args[@]}"; else python -u -m "$1" "${args[@]}"; fi
echo "end $(date)"
