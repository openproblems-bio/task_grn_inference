#!/bin/bash
#SBATCH --job-name=workflow_local
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=10:00:00
#SBATCH --mem=250GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com

set -e

DATASETS=(
    "op"
)
METHODS=(
    "scprint"
)

# - where to save the scores (all metrics, datasets, methods)
SAVE_SCORES_FILE="output/scores.csv"
# - whether to force re-run the inference if the files exists
FORCE=true
# - whether to run the consensus for reg2 (only run when to update the consensus)
RUN_CONSENSUS_FLAG=False

# ----- run methods -----
cmd="python src/workflows_local/benchmark/methods/script.py 
        --datasets ${DATASETS[@]} 
        --methods ${METHODS[@]}"

[ "$FORCE" = true ] && cmd="${cmd} --force"

echo "Running: $cmd"
$cmd

# # ----- run metrics -----
# cmd="python src/workflows_local/benchmark/metrics/script.py 
#         --datasets ${DATASETS[@]} 
#         --methods ${METHODS[@]}
#         --save_scores_file ${SAVE_SCORES_FILE}"

# [ "$RUN_CONSENSUS_FLAG" = true ] && cmd="${cmd} --run_consensus_flag"

# echo "Running: $cmd"
# $cmd