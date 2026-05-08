
#!/bin/bash
set -e

run_prefix=${run_prefix:-'sbatch'}
python src/utils/config.py
source src/utils/config.env

# DATASETS=(${DATASETS//,/ })
DATASETS=('soundlife' 'soundlife_vaccine') #'op' 'adamson' 'replogle' 'norman' 'nakatake' 'parsebioscience'  '300BCG' 'xaira_HCT116' 'xaira_HEK293T' 'ibd_uc' 'ibd_cd'

# METHODS=(${METHODS//,/ })
METHODS=('geneformer')

methods_dir='src/methods/'
ctr_methods_dir='src/methods/'
resources_dir='resources'
predictions_dir="$resources_dir/results/benchmark/"


MULTIOMICS_METHODS=("scenicplus" "scglue" "figr" "granie" "celloracle")

run_func() {
    local method=$1
    local dataset=$2

    rna="${resources_dir}/grn_benchmark/inference_data/${dataset}_rna.h5ad"
    rna_all="${resources_dir}/extended_data/${dataset}_rna_all.h5ad"
    atac="${resources_dir}/grn_benchmark/inference_data/${dataset}_atac.h5ad"
    prediction="$predictions_dir/${dataset}/${dataset}.${method}.${method}.prediction.h5ad"

    if [[ ! -f "$atac" ]] && [[ " ${MULTIOMICS_METHODS[*]} " =~ " ${method} " ]]; then
        echo "Skipping $method on $dataset (no ATAC data)"
        return
    fi

    echo "Running $method on $dataset"
    script="${methods_dir}${method}/run_local.sh"
    arguments="--rna $rna --prediction $prediction --atac $atac --rna_all $rna_all"

    if [[ "$run_prefix" == "bash" ]]; then
        bash "$script" $arguments
    elif [[ "$run_prefix" == "sbatch" ]]; then
        output=$(sbatch "$script" $arguments)
        echo "$output"
        jobid=$(echo "$output" | awk '{print $4}')
        echo "Job ID: $jobid"
    else
        echo "Unknown run_prefix: $run_prefix"
    fi
}

for dataset in "${DATASETS[@]}"; do
    mkdir "$predictions_dir/$dataset" || true
    for method in "${METHODS[@]}"; do
        run_func "$method" "$dataset"
    done
done