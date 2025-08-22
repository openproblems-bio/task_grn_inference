set -e

datasets=('xaira_HEK293T') #'replogle' 'op' 'nakatake' 'adamson' 'norman'  'xaira_HEK293T' 'parsescience'
run_local=true # set to true to run locally, false to run on AWS

run_grn_inference=true
run_grn_evaluation=false
run_download=false


for dataset in "${datasets[@]}"; do

    if [ "$run_grn_inference" = true ]; then
        echo "Running GRN inference for dataset: $dataset"
        if [ "$run_local" = true ]; then
            echo "Running locally"
            
            file="resources/results/$dataset/trace.txt"

            if [ -f "$file" ]; then
                
                dir=$(dirname "$file")
                base=$(basename "$file" .txt)
                today=$(date +%Y-%m-%d)
                cp "$file" "${dir}/${base}_${today}.txt"
            fi
        else
            echo "Running on AWS"
        fi
        bash scripts/run_grn_inference.sh --dataset=$dataset --run_local=$run_local
        
    fi

    if [ "$run_grn_evaluation" = true ]; then
        if [ "$run_local" = false ]; then
            
            file="resources/results/$dataset/trace.txt"

            if [ -f "$file" ]; then
                echo "Making a copy of previous trace file"
                dir=$(dirname "$file")
                base=$(basename "$file" .txt)
                today=$(date +%Y-%m-%d)
                cp "$file" "${dir}/${base}_${today}.txt"
            fi

            echo "Downloading inference results from AWS"
            aws s3 sync  s3://openproblems-data/resources/grn/results/$dataset resources/results/$dataset 
        fi 
        echo "Running consensus for dataset: $dataset"
        bash scripts/prior/run_consensus.sh $dataset # run consensus for regression 2 and ws distance -> needs to be run after adding each method and dataset
        
        if [ "$run_local" = false ]; then
            echo "Syncing prior results to AWS"
            aws s3 sync  resources/grn_benchmark/prior s3://openproblems-data/resources/grn/prior 
        fi

        echo "Running GRN evaluation for dataset: $dataset"
        bash scripts/run_grn_evaluation.sh --dataset=$dataset --run_local=$run_local --build_images=false 
    fi

    if [ "$run_download" = true ]; then
        if [ "$run_local" = false ]; then
            echo "Downloading evaluation results from AWS"
            aws s3 sync  s3://openproblems-data/resources/grn/results/$dataset resources/results/$dataset 
        fi
    fi
    # bash scripts/prior/run_ws_background.sh # run background distance for ws distance -> needs to be run after adding each dataset

done




