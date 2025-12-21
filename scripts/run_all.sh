set -e

python src/utils/config.py
source src/utils/config.env
DATASETS=(${DATASETS//,/ })

run_local=false
run_grn_inference=false #arg
run_consensus=true
run_grn_evaluation=true #arg
run_sync=false

num_workers=20


for dataset in "${DATASETS[@]}"; do
    trace_file="resources/results/$dataset/trace.txt"

    if [ "$run_grn_inference" = true ]; then
        echo "Running GRN inference for dataset: $dataset"
        if [ "$run_local" = true ]; then
            echo "Running locally"

            if [ -f "$trace_file" ]; then
                dir=$(dirname "$trace_file")
                base=$(basename "$trace_file" .txt)
                today=$(date +%Y-%m-%d)
                cp "$trace_file" "${dir}/${base}_${today}.txt"
            fi
        else
            echo "Running on AWS"
        fi
        bash scripts/run_grn_inference.sh --dataset=$dataset --run_local=$run_local
        
    fi

    if [ "$run_grn_evaluation" = true ]; then

        if [ -f "$trace_file" ]; then
            dir=$(dirname "$trace_file")
            base=$(basename "$trace_file" .txt)
            today=$(date +%Y-%m-%d)
            cp "$trace_file" "${dir}/${base}_${today}.txt"
        fi
        
        # if [ "$run_local" = false ]; then
        #     echo "Uploading inference results to AWS"
        #     aws s3 sync  resources/results/$dataset s3://openproblems-data/resources/grn/results/$dataset 
        #     aws s3 sync  s3://openproblems-data/resources/grn/results/$dataset resources/results/$dataset 
        # fi 
        

        # if [ "$run_local" = false ]; then
        #     echo "Downloading inference results from AWS"
        #     aws s3 sync  s3://openproblems-data/resources/grn/results/$dataset resources/results/$dataset 
        # fi
        if [ "$run_consensus" = true ]; then
            echo "Running consensus for dataset: $dataset"
            bash scripts/prior/run_consensus.sh --dataset $dataset # run consensus for Regression and ws distance -> needs to be run after adding each method and dataset
        fi
        
        if [ "$run_local" = false ]; then
            echo "Syncing prior results to AWS"
            aws s3 sync  resources/grn_benchmark/prior s3://openproblems-data/resources/grn/grn_benchmark/prior 
        fi

        echo "Running GRN evaluation for dataset: $dataset"
        bash scripts/run_grn_evaluation.sh --dataset=$dataset --run_local=$run_local --build_images=false --num_workers=$num_workers
    fi

    if [ "$run_sync" = true ]; then
        if [ "$run_local" = false ]; then
            echo "Downloading evaluation results from AWS"
            aws s3 sync  s3://openproblems-data/resources/grn/results/$dataset resources/results/$dataset 
        fi
    fi

    # bash scripts/prior/run_ws_background.sh # run background distance for ws distance -> needs to be run after adding each dataset

done