
import os
import sys
import subprocess
import anndata as ad
import numpy as np


def run_grn_inference(dataset='op', subsample=None):
    par = {
        'models_dir': f'resources/grn_models/{dataset}/',
        'rna': f'resources/inference_datasets/{dataset}_rna.h5ad',
        'atac': f'resources/inference_datasets/{dataset}_atac.h5ad', 
        'rna_positive_control': f'resources/datasets_raw/{dataset}.h5ad', 
        'num_workers': 10,
        'tmp_dir': 'output/grn_inference'
        # 'max_n_links': 10000,
    }
   
    os.makedirs(par['models_dir'], exist_ok=True)
    os.makedirs(par['tmp_dir'], exist_ok=True)

    if subsample is not None:
        if dataset=='op':
            new_atac = f"{par['tmp_dir']}/{dataset}_atac.h5ad"
            new_rna = f"{par['tmp_dir']}/{dataset}_rna.h5ad"
            if not os.path.exists(new_rna):
                rna = ad.read_h5ad(par['rna'])
                atac = ad.read_h5ad(par['atac'])
                if subsample == 1:
                    mask_rna = rna.obs['donor_id'] == 'donor_0'
                    mask_atac = atac.obs['donor_id'] == 'donor_0'
                elif subsample == 2:
                    mask_rna = rna.obs['donor_id'].isin(['donor_0', 'donor_1'])
                    mask_atac = atac.obs['donor_id'].isin(['donor_0', 'donor_1'])
                rna = rna[mask_rna, :]
                atac = atac[mask_rna, :]
                
                rna.write(new_rna)
                atac.write(new_atac)
                print(rna)
            par['rna'] = new_rna
            par['atac'] = new_atac
        else:
            new_rna = f"{par['tmp_dir']}/{dataset}_rna.h5ad"
            if not os.path.exists(new_rna):
                rna = ad.read_h5ad(par['rna'])
                all_perturbation = rna.obs['perturbation'].unique()
                sample_size = int(len(all_perturbation) * subsample)  # Number of samples
                selected_samples = np.random.choice(all_perturbation, sample_size, replace=False)
                rna = rna[rna.obs['perturbation'].isin(selected_samples), :]
                rna.write(new_rna)
                print(rna)
            par['rna'] = new_rna
        

    for method in methods:
        print(method)
        if subsample is None:
            par['prediction'] = f"{par['models_dir']}/{method}.csv"
        else:
            par['prediction'] = f"{par['models_dir']}/{method}_{subsample}.csv"
        
        if (force == False) & (os.path.exists(par['prediction'])):
            continue
        
        # Method arguments
        if method == 'positive_control':
            method_args = (f"--rna {par['rna_positive_control']} ")
        else:
            method_args = (f"--rna {par['rna']} ")

        method_args += (
                        f"--prediction {par['prediction']} "
                        f"--num_workers {par['num_workers']} "
                        #    f"--max_n_links {par['max_n_links']} "
                        )
        
        method_args_r = ( f" {par['prediction']} "
                        
        )
        # Determine the command based on the method
        if method in ["positive_control", "pearson_corr", "negative_control"]:
            command = f"python src/control_methods/{method}/script.py {method_args}"
        elif method == "celloracle":
            method_args += f"--atac {par['atac']} "
            command = (f"/home/jnourisa/miniconda3/envs/celloracle/bin/python "
                       f"src/methods/multi_omics/celloracle/script.py {method_args}")
        elif method in ["grnboost2", "scenic", "genie3"]:
            command = f"singularity exec ../../images/scenic python src/methods/single_omics/{method}/script.py {method_args}"
        elif method == 'scglue':
            method_args += f"--atac {par['atac']} "
            command = f"singularity exec ../../images/scglue python src/methods/multi_omics/{method}/script.py {method_args}"
        elif method == 'scenicplus':
            method_args += f"--atac {par['atac']} "
            command = f"singularity exec ../../images/scenicplus python src/methods/multi_omics/{method}/script.py {method_args}"
        elif method=='ppcor':
            command = f"singularity exec ../../images/ppcor Rscript src/methods/single_omics/{method}/script.R {method_args}"
        else:
            command = f"singularity exec ../../images/{method} python src/methods/single_omics/{method}/script.py {method_args}"
        

        if method in ["positive_control", "pearson_corr", "negative_control"]:
            mem = "64GB"
            time = "0:30:00"
        elif method in ["grnboost2"]:
            mem = "120GB"
            time = "15:00:00"
        elif method in ["portia"]:
            mem = "250GB"
            time = "6:00:00"
        elif method in ["celloracle"]:
            mem = "250GB"
            time = "12:00:00"
        elif method in ["scenicplus"]:
            mem = "250GB"
            time = "12:00:00"
        elif method in ["scenic"]:
            mem = "250GB"
            time = "24:00:00"

        # Prepare sbatch command
        tag = f"--job-name={method}"  # No spaces around '='
        resources = (f"--cpus-per-task={par['num_workers']} "
                     f"--mem={mem} --time={time} --partition={partition}")
        
        # Combine tags and resources
        full_tag = [tag] + resources.split()
        
        # Add GPU partition if method is 'scglue'
        if method == 'scglue':
            full_tag += ["--partition=gpu", "--gres=gpu:1"]

        try:
            if sbatch:
                result = subprocess.run(['sbatch'] + full_tag + ['scripts/sbatch/grn_inference.sh', command], check=True, capture_output=True, text=True)
            else:
                result = subprocess.run(['bash'] + ['scripts/sbatch/grn_inference.sh', command], check=True, capture_output=True, text=True)

            print(f"Job {method} submitted successfully.")
            print(result.stdout)  # Print the standard output
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while submitting job for {method}: {e}")
            print(f"Return code: {e.returncode}")
            print(f"Command output: {e.output}")
            print(f"Command error output: {e.stderr}")

if __name__ == '__main__':
    force = True
    sbatch = True
    # methods = ["positive_control", "negative_control", "pearson_corr", "portia", "grnboost2", "ppcor", "scenic"],
    # methods = ["portia", "grnboost2"]
    methods = ["scenic"]
    datasets = ['adamson']

    
    partition='cpu'

    # mem = "120GB"
    # time = "24:00:00"

    if True: # normal run 
        for dataset in datasets:
            run_grn_inference(dataset, subsample=None)

    if False: # subsample
        # for dataset in ['replogle2', 'norman', 'adamson', 'nakatake']: # 'replogle2' 'op' norman
        for dataset in ['op']:
            if dataset == 'op':
                for subsample in [1, 2]: #number of donors 
                    run_grn_inference(dataset, subsample=subsample)
            else:
                for subsample in [0.2, 0.5, 1.0]:
                    run_grn_inference(dataset, subsample=subsample)