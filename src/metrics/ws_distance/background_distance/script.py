import anndata as ad
from tqdm import tqdm
import scipy
import pandas as pd
import numpy as np
import multiprocessing as mp
import os

# Shared memory path for NumPy array
shared_mem_path = "/home/jnourisa/projs/ongoing/ciim/output/adata_X.npy"

def load_adata(par):
    """Load AnnData and store X in a shared memory-mapped NumPy array."""
    adata = ad.read_h5ad(par['evaluation_data_sc'])
    if 'X_norm' in adata.layers:
        layer = 'X_norm'
    elif 'lognorm' in adata.layers:
        layer = 'lognorm'
    elif 'pearson_residual' in adata.layers:
        layer = 'pearson_residual'
    else:
        raise ValueError("No suitable layer found in AnnData. Please provide a valid layer.")
    adata.X = adata.layers[layer]

    if 'is_control' not in adata.obs.columns:
        adata.obs['is_control'] = adata.obs['perturbation'] == 'non-targeting'

    # Convert to a shared memory-mapped NumPy array
    np.save(shared_mem_path, adata.X.toarray())  # Save sparse matrix as dense NumPy
    return adata

def calculate_ws_distance(net, gene_names, obs_perturbation, progress):
    """
    Compute Wasserstein distance for gene pairs in `net`.
    Each worker loads the shared memory-mapped array instead of the full `adata`.
    """
    ws_distances = []
    
    # Load shared memory-mapped NumPy array
    X_shared = np.load(shared_mem_path, mmap_mode='r')

    for parent, child in zip(net['source'], net['target']):
        mask_gene = np.where(gene_names == child)[0]

        mask_sample_ctrl = obs_perturbation == 'control'
        mask_sample_intv = obs_perturbation == parent

        X_observ = X_shared[np.ix_(mask_sample_ctrl, mask_gene)]
        X_interv = X_shared[np.ix_(mask_sample_intv, mask_gene)]

        assert X_observ.shape[0] != 0
        assert X_interv.shape[0] != 0

        ws_distance = scipy.stats.wasserstein_distance(
            X_observ.reshape(-1), X_interv.reshape(-1)
        )

        ws_distances.append(ws_distance)

    net['ws_distance'] = ws_distances

    # Update progress counter directly
    progress.value += 1

    return net

def process_tf(tf, gene_names, obs_perturbation, progress):
    """Compute Wasserstein distance for all genes targeted by a TF."""
    net_all = pd.DataFrame({'source': tf, 'target': gene_names})
    return calculate_ws_distance(net_all, gene_names, obs_perturbation, progress)

import argparse
arg = argparse.ArgumentParser(description='Compute Wasserstein distance for background GRN inference')
arg.add_argument('--dataset', type=str, help='Dataset to use for the analysis')
arg.add_argument('--evaluation_data_sc', type=str, help='Path to the evaluation data in single-cell format')
arg.add_argument('--background_distance', type=str, help='Path to save the background distance results')
arg.add_argument('--tf_all', type=str, help='Path to the file containing all transcription factors')
arg.add_argument('--layer', type=str, default='X_norm', help='Layer to use for the analysis')
arg.add_argument('--num_workers', type=int, default=100, help='Maximum number of workers for multiprocessing')
args = arg.parse_args()
par = args.__dict__


def main(par):
    
    # Load AnnData once and store in shared memory
    adata = load_adata(par)

    # Extract metadata to pass to workers
    adata_shape = adata.shape
    gene_names = adata.var_names.to_numpy()
    
    adata.obs['perturbation'] = adata.obs['perturbation'].astype(str)
    if 'is_control' not in adata.obs.columns:
        adata.obs['is_control'] = adata.obs['perturbation'] == 'non-targeting'
    adata.obs.loc[adata.obs['is_control'], 'perturbation'] = 'control'

    obs_perturbation = adata.obs['perturbation'].to_numpy()

    # Get available TFs
    tf_all = np.loadtxt(par['tf_all'], dtype=str)
    available_tfs = np.intersect1d(obs_perturbation, tf_all)

    num_workers = min(par['num_workers'], len(available_tfs))

    # Use multiprocessing with a progress tracker
    with mp.Manager() as manager:
        progress = manager.Value('i', 0)  # Shared integer counter

        with mp.Pool(num_workers) as pool:
            results = []
            with tqdm(total=len(available_tfs), desc="Processing TFs", position=0) as pbar:
                for result in pool.starmap_async(
                    process_tf,
                    [(tf, gene_names, obs_perturbation, progress) for tf in available_tfs]
                ).get():
                    results.append(result)
                    pbar.update(1)  # Update tqdm after each job completes
                    pbar.refresh()

    # Save results
    net_all_ws_distance = pd.concat(results, ignore_index=True)
    net_all_ws_distance.to_csv(par['background_distance'], index=False)

    # Clean up shared memory
    os.remove(shared_mem_path)

if __name__ == '__main__':
    main(par)