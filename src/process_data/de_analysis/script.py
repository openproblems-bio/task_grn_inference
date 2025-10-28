import scanpy as sc
import pandas as pd
import numpy as np
import scanpy as sc
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm
import os
import sys


try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
    "util_dir": 'src/utils'
    }
    sys.path.append(meta["util_dir"])
from util import manage_layer


for dataset in ["replogle", "xaira_HEK293T", "xaira_HCT116"]:
    single_cell_file = f'resources/processed_data/{dataset}_evaluation_sc.h5ad'

    if not os.path.exists(single_cell_file):
        print(f"Skipping {dataset}, file not found.")
        continue
    else:
        print(f'------ {dataset} ------ ')

    # load data
    adata = sc.read_h5ad(single_cell_file, backed='r')
    perturbation_type = adata.uns['perturbation_type']
    tf_all = pd.read_csv('resources/grn_benchmark/prior/tf_all.csv', header=None)[0].tolist()

    # subset to only TF perturbations
    pert_mask = adata.obs['perturbation'].isin(tf_all)
    ctrl_mask = adata.obs['is_control']

    # subsample controls to 10,000 if needed
    ctrl_indices = np.where(ctrl_mask)[0]
    if len(ctrl_indices) > 10000:
        np.random.seed(0)  # for reproducibility
        ctrl_indices = np.random.choice(ctrl_indices, size=10000, replace=False)

    # combine perturbations + subsampled controls
    subset_indices = np.concatenate([np.where(pert_mask)[0], ctrl_indices])
    adata = adata[subset_indices, :].to_memory()
    assert adata.shape[0] > 0

    layer = manage_layer(adata, par={'layer': 'lognorm'}) 

    adata.X = adata.layers[layer]
    adata.X = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X



    def wrapper_de_analysis(adata):
        # Cache control subset
        ctrl = adata[adata.obs['is_control']].copy()
        assert ctrl.n_obs > 0, "No control cells found!"

        # Prepare mapping from TF -> indices in adata
        idx_dict = {}
        for tf in tf_all:
            pert_mask = adata.obs['perturbation'] == tf
            if pert_mask.sum() == 0:
                continue
            idx_dict[tf] = {
                'pert_idx': np.where(pert_mask)[0],
                'ctrl_idx': np.where(adata.obs['is_control'])[0]
            }

        def compute_de_for_tf_indices(tf, idx_dict, X, obs_index, var):
            pert_idx = idx_dict['pert_idx']
            ctrl_idx = idx_dict['ctrl_idx']

            adata_sub = sc.AnnData(
                X=np.vstack([X[pert_idx, :], X[ctrl_idx, :]]),
                obs=pd.concat([
                    pd.DataFrame({"group":"pert"}, index=obs_index[pert_idx]),
                    pd.DataFrame({"group":"ctrl"}, index=obs_index[ctrl_idx])
                ]),
                var=var.copy()
            )
            adata_sub.obs['group'] = adata_sub.obs['group'].astype('category')

            sc.tl.rank_genes_groups(
                adata_sub,
                groupby='group',
                groups=['pert'],
                reference='ctrl',
                method='wilcoxon',
                # n_genes=adata_sub.n_vars,
                use_raw=False
            )
            df = sc.get.rank_genes_groups_df(adata_sub, group='pert')
            # create wide format: one row per perturbation (TF), columns are genes
            df_wide = df.set_index('names')['scores'].to_frame().T  # make it a single-row df
            df_wide.index = [tf]  # set the TF as the row index
            df_wide['perturbation'] = tf
            df_wide = df_wide.set_index('perturbation')
            return df_wide

        # Run in parallel
        de_results = Parallel(n_jobs=10, backend='loky', verbose=10)(
            delayed(compute_de_for_tf_indices)(tf, idx_dict, adata.X, adata.obs.index, adata.var[[]])
            for tf, idx_dict in tqdm(idx_dict.items(), desc="Processing TFs")
        )

        de_all_wide = pd.concat([r for r in de_results if r is not None], axis=0)

        return de_all_wide
    de_all = wrapper_de_analysis(adata)

    X_scores = de_all.values  # shape: n_TFs × n_genes
    obs = pd.DataFrame(index=de_all.index)  # each row is a TF
    var = pd.DataFrame(index=de_all.columns)  # each column is a gene

    # create AnnData
    adata_de = sc.AnnData(X=X_scores, obs=obs, var=var)

    print(adata_de)

    adata_de.uns['dataset_id'] = dataset
    adata_de.uns['perturbation_type'] = perturbation_type

    for ky, value in adata.uns.items():
        adata_de.uns[ky] = value

    adata_de.write(f'resources/grn_benchmark/evaluation_data/{dataset}_de.h5ad')
