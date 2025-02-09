import anndata as ad
from tqdm import tqdm
import scipy
import pandas as pd
import numpy as np
import multiprocessing as mp
import os

# Configuration parameters
par = {
    'evaluation_data_sc': 'resources/grn_benchmark/evaluation_data/replogle_sc.h5ad',
    'background_distance': 'resources/grn_benchmark/prior/ws_distance_background_replogle.csv',
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
    'layer': 'X_norm',
    'max_workers': 100
}

# Shared memory path for NumPy array
shared_mem_path = "/tmp/adata_X.npy"

def load_adata(par):
    """Load AnnData and store X in a shared memory-mapped NumPy array."""
    adata = ad.read_h5ad(par['evaluation_data_sc'])
    adata.X = adata.layers[par['layer']]

    if 'is_control' not in adata.obs.columns:
        adata.obs['is_control'] = adata.obs['perturbation'] == 'non-targeting'

    # Convert to a shared memory-mapped NumPy array
    np.save(shared_mem_path, adata.X.toarray())  # Save sparse matrix as dense NumPy
    return adata

def calculate_ws_distance(net, adata_shape, gene_names, obs_perturbation, progress):
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

def process_tf(tf, adata_shape, gene_names, obs_perturbation, progress):
    """Compute Wasserstein distance for all genes targeted by a TF."""
    net_all = pd.DataFrame({'source': tf, 'target': gene_names})
    return calculate_ws_distance(net_all, adata_shape, gene_names, obs_perturbation, progress)

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

    num_workers = min(par['max_workers'], len(available_tfs))

    # Use multiprocessing with a progress tracker
    with mp.Manager() as manager:
        progress = manager.Value('i', 0)  # Shared integer counter

        with mp.Pool(num_workers) as pool:
            results = []
            with tqdm(total=len(available_tfs), desc="Processing TFs", position=0) as pbar:
                for result in pool.starmap_async(
                    process_tf,
                    [(tf, adata_shape, gene_names, obs_perturbation, progress) for tf in available_tfs]
                ).get():
                    results.append(result)
                    pbar.update(1)  # Update tqdm after each job completes

    # Save results
    net_all_ws_distance = pd.concat(results, ignore_index=True)
    net_all_ws_distance.to_csv(par['background_distance'], index=False)

    # Clean up shared memory
    os.remove(shared_mem_path)

if __name__ == '__main__':
    main(par)