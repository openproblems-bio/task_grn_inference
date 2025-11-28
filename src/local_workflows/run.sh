
#!/bin/bash
set -e

run_prefix='sbatch' #bash
DATASETS=('op' 'adamson' 'replogle' 'norman' 'nakatake' 'parsebioscience'  '300BCG' 'xaira_HCT116' 'xaira_HEK293T') #'op' 'adamson' 'replogle' 'norman' 'nakatake' 'parsebioscience'  '300BCG' 'xaira_HCT116' 'xaira_HEK293T' 'ibd_uc' 'ibd_cd'
DATASETS=('parsebioscience') #'op' 'adamson' 'replogle' 'norman' 'nakatake' 'parsebioscience'  '300BCG' 'xaira_HCT116' 'xaira_HEK293T' 'ibd_uc' 'ibd_cd'

# METHODS=('negative_control' 'positive_control' 'pearson_corr' 'portia' 'ppcor' 'grnboost' 'scenic'  'scenicplus' 'scglue' 'figr' 'granie')
METHODS=( 'scenic' 'grnboost') #'negative_control' 'positive_control' 'pearson_corr' 'portia' 'ppcor' 'grnboost' 'scenic'  'scenicplus' 'scglue' 'figr' 'granie'

methods_dir='src/methods/'
ctr_methods_dir='src/methods/'
resources_dir='resources'
predictions_dir="$resources_dir/results/"


run_func() {
    local method=$1
    local dataset=$2

    echo "Running $method on $dataset"

    script="${methods_dir}${method}/run_local.sh"

    rna="${resources_dir}/grn_benchmark/inference_data/${dataset}_rna.h5ad"
    rna_all="${resources_dir}/extended_data/${dataset}_bulk.h5ad"
    atac="${resources_dir}/grn_benchmark/inference_data/${dataset}_atac.h5ad"
    prediction="$predictions_dir/${dataset}/${dataset}.${method}.${method}.prediction.h5ad"

    arguments="--rna $rna --prediction $prediction --atac $atac --rna_all $rna_all"

    if [[ "$run_prefix" == "bash" ]]; then
        bash "$script" $arguments
    elif [[ "$run_prefix" == "sbatch" ]]; then
        # submit the job and capture the job ID
        output=$(sbatch "$script" $arguments)
        echo "$output"
        # sbatch usually returns: "Submitted batch job 12345678"
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