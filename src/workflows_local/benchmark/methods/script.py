
import os
import sys
import subprocess
import anndata as ad
import numpy as np
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument('--force', action='store_true', help='Force overwrite exisiting GRN predictions')
argparser.add_argument('--sbatch', action='store_true', help='Submit jobs to SLURM')
argparser.add_argument('--partition', type=str, default='cpu', help='Partition to submit jobs to')
argparser.add_argument('--methods', type=str, nargs='+', required=True,help='Methods to run')
argparser.add_argument('--datasets', type=str, nargs='+', required=True, help='Datasets to run inference on')
args = argparser.parse_args()
par = vars(args)



def get_par(dataset, method):
    par_local = {
        'models_dir': f'resources/grn_models/{dataset}/',
        'rna': f'resources/grn_benchmark/inference_data/{dataset}_rna.h5ad',
        'atac': f'resources/grn_benchmark/inference_data/{dataset}_atac.h5ad', 
        'dataset_id': dataset,
        'num_workers': 10,
        'tmp_dir': 'output/grn_inference'
    }

    for key, value in par_local.items():
        par[key] = value

    if (method == 'scprint') & (dataset == 'replogle'):
        par['rna'] = f'resources/grn_benchmark/inference_data/replogle_rna_sc_subset.h5ad'
    if (method == 'positive_control') & (dataset == 'op'):
        par['rna'] = f'resources/grn_benchmark/evaluation_data/{dataset}_bulk.h5ad'
    if (method == 'positive_control') & (dataset in ['replogle', 'norman', 'adamson', 'nakatake']):
        par['rna'] = f'resources/extended_data/{dataset}_bulk.h5ad'

    return par

def run_grn_inference(par, dataset='op', subsample=None):
    
    for method in par['methods']:
        
        if (method=='scprint') & (dataset in ['nakatake', 'adamson']):
            continue
        
        par = get_par(dataset, method)

        os.makedirs(par['models_dir'], exist_ok=True)
        os.makedirs(par['tmp_dir'], exist_ok=True)

        if subsample is None:
            par['prediction'] = f"{par['models_dir']}/{method}.h5ad"
        else:
            par['prediction'] = f"{par['models_dir']}/{method}_{subsample}.h5ad"
        
        if (par['force'] == False) & (os.path.exists(par['prediction'])):
            continue
        
        # Method arguments

        method_args = (
                        f"--rna {par['rna']} "
                        f"--prediction {par['prediction']} "
                        f"--num_workers {par['num_workers']} "
                        f"--dataset_id {par['dataset_id']} "
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
        elif method=='scprint':
            # command = "source ~/miniconda3/bin/activate scprint \n"
            command = f"python src/methods/single_omics/{method}/script.py {method_args}"
        elif method=='figr':
            # command = "source ~/miniconda3/bin/activate scprint \n"
            command = f"singularity exec ../../images/figr Rscript src/methods/multi_omics/figr/script.R {method_args}"
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
        elif method in ["scprint"]:
            if dataset == 'replogle':
                mem = "500GB"
                time = "24:00:00"
            else:
                mem = "120GB"
                time = "24:00:00"
        elif method in ["scglue", "figr"]:
            mem = "250GB"
            time = "24:00:00"


        else:
            # print(f"Method {method} not recognized")
            raise ValueError(f"Method {method} not recognized")

        # Prepare sbatch command
        tag = f"--job-name={method}"  # No spaces around '='
        resources = (f"--cpus-per-task={par['num_workers']} "
                     f"--mem={mem} --time={time} --partition={par['partition']}")
        
        # Combine tags and resources
        full_tag = [tag] + resources.split()
        
        # Add GPU partition if method is 'scglue'
        if method in ['scglue', 'scprint']:
            full_tag += ["--partition=gpu", "--gres=gpu:1"]

        try:
            if par['sbatch']:
                result = subprocess.run(['sbatch'] + full_tag + ['scripts/sbatch/grn_inference.sh', command], check=True, capture_output=True, text=True)
            else:
                result = subprocess.run(['bash'] + ['scripts/sbatch/grn_inference.sh', command], check=True, capture_output=True, text=True)

            if result.returncode != 0:
                print("STDOUT:", result.stdout)
                print("STDERR:", result.stderr)
                raise RuntimeError(f"Error: command process dataset failed with exit code {result.returncode}")
            print(result.stdout)  # Print the standard output
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while submitting job for {method}: {e}")
            print(f"Return code: {e.returncode}")
            print(f"Command output: {e.output}")
            print(f"Command error output: {e.stderr}")

def main(par):
    for dataset in par['datasets']:
        run_grn_inference(par, dataset, subsample=None)


if __name__ == '__main__':
    main(par)
    
    # if False: # subsample
    #     # for dataset in ['replogle', 'norman', 'adamson', 'nakatake']: # 'replogle' 'op' norman
    #     for dataset in ['op']:
    #         if dataset == 'op':
    #             for subsample in [1, 2]: #number of donors 
    #                 run_grn_inference(dataset, subsample=subsample)
    #         else:
    #             for subsample in [0.2, 0.5, 1.0]:
    #                 run_grn_inference(dataset, subsample=subsample)