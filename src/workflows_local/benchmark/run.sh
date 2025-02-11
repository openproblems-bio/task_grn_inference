#!/bin/bash
#SBATCH --job-name=workflow_local
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00
#SBATCH --mem=500GB
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com

set -e
source ~/miniconda3/bin/activate scprint

# ----- parameters -----
RUN_GRN_INFERENCE=true
CALCULATE_METRICS=true
RUN_CONSENSUS_FLAG=false # - whether to run the consensus for reg2 (only run when to update the consensus), #TODO: update the code to handle other consensus
FORCE=true
SBATCH=false
SAVE_SCORES_FILE="resources/scores/scores_controls.csv" # - where to save the scores (all metrics, datasets, methods)


DATASETS=(
    " replogle op norman nakatake adamson "
)

# DATASETS=(
#     " op "
# )

# METHODS=(
#     "scglue scenicplus celloracle granie figr collectri grnboost2 ppcor portia scenic positive_control pearson_corr negative_control scprint"
# )

METHODS=(
    " pearson_corr positive_control "
)


if [ "$RUN_GRN_INFERENCE" = true ]; then
    # ----- run grn inference -----
    cmd="python src/workflows_local/benchmark/methods/script.py \
        --datasets ${DATASETS[@]} \
        --methods ${METHODS[@]}"

    [ "$FORCE" = true ] && cmd="${cmd} --force"
    [ "$SBATCH" = true ] && cmd="${cmd} --sbatch"

    echo "Running: $cmd"
    eval "$cmd"
fi
if [ "$CALCULATE_METRICS" = true ]; then
    # ----- run metrics -----
    if [ "$SBATCH" = false ]; then
        # ----- run metrics -----
        cmd="python src/workflows_local/benchmark/metrics/script.py 
                --datasets ${DATASETS[@]} 
                --methods ${METHODS[@]}
                --save_scores_file ${SAVE_SCORES_FILE}"

        [ "$RUN_CONSENSUS_FLAG" = true ] && cmd="${cmd} --run_consensus_flag"

        echo "Running: $cmd"
        $cmd
    fi

    echo $SAVE_SCORES_FILE
fi


