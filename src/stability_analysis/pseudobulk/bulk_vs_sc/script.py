import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import sys

meta = {
    "resources_dir": './',
    "util_dir": 'src/utils'
}
sys.path.append(meta["resources_dir"])
sys.path.append(meta["util_dir"])

# from src.control_methods.pearson_corr.script import main as main_inference
from src.metrics.regression_2.helper import main as main_reg2
from src.metrics.ws_distance.helper import main as main_ws_distance
from src.metrics.ws_distance.consensus.helper import main as main_consensus_ws_distance


def def_par(dataset):
    par = {
        'evaluation_data': f'resources/grn_benchmark/evaluation_data/{dataset}_bulk.h5ad',
        'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
        'apply_skeleton': False,
        'apply_tf': True,
        'max_n_links': 50000,
        'layer': 'lognorm',
        'apply_tf_methods': True,
        'reg_type': 'ridge',
        'num_workers': 10,
        'ws_consensus': f'resources/grn_benchmark/prior/ws_consensus_{dataset}.csv',
        'ws_distance_background': f'resources/grn_benchmark/prior/ws_distance_background_{dataset}.csv',
        'evaluation_data_sc': f'resources/processed_data/{dataset}_sc.h5ad'

    }
    return par
def prediction_file_name(dataset, data_type):
    return f'{results_dir}/{dataset}.prediction_{data_type}.h5ad'

def infer_grn(par, dataset):
    from util import corr_net
    adata = ad.read_h5ad(par["rna"], backed='r')
    tf_all = np.loadtxt(par["tf_all"], dtype=str)

    if dataset == 'parsebioscience':
        perturbs = sorted(adata.obs['perturbation'].unique())[:10]
    else:
        perturbs = tf_all
    adata = adata[adata.obs['perturbation'].isin(perturbs)]
    print(adata.shape)
    adata = adata.to_memory()

    net = corr_net(adata, tf_all, par)
    net = net.astype(str)
    net = ad.AnnData(
        X=None,
        uns={
            "method_id": 'pearson_corr',
            "dataset_id": adata.uns['dataset_id'],
            "prediction": net[["source", "target", "weight"]]
        }
    )
    return net

if __name__ == '__main__':
    results_dir = 'resources/results/experiment/bulk_vs_sc'
    os.makedirs(results_dir, exist_ok=True)
    datasets = ['parsebioscience'] #'replogle', 'xaira_HEK293T', 'parsebioscience' 'xaira_HCT116'

    metrics_all = []
    for dataset in datasets:
        par = def_par(dataset)
        print('Processing dataset:', dataset, flush=True)
        if True:
            # - infer GRNs
            for data_type in ['sc', 'bulk']: 
                print(f"Inferring GRNs for {data_type} data...", flush=True)
                if data_type == 'bulk':
                    par['rna'] = f'resources/grn_benchmark/inference_data/{dataset}_rna.h5ad'
                else:
                    par['rna'] = f'resources/extended_data/{dataset}_train_sc.h5ad'
                net = infer_grn(par, dataset)

                par['prediction'] = prediction_file_name(dataset, data_type)
                net.write_h5ad(par['prediction'])
        
        # - consensus 
        par['regulators_consensus'] = f'{results_dir}/regulators_consensus_{dataset}.json'
        par['ws_consensus'] = f'{results_dir}/ws_consensus_{dataset}.json'
        if True:
            def naming_convention(dataset, model):
                return f'{dataset}.{model}.h5ad'
            from src.metrics.regression_2.consensus.helper import main as main_consensus_reg2
            par['naming_convention'] = naming_convention
            par['dataset'] = dataset
            par['models_dir'] = results_dir
            par['models'] = ['prediction_bulk', 'prediction_sc']
            _ = main_consensus_reg2(par)
            main_consensus_ws_distance(par)

        # - grn evaluation
        rr_all_store = []
        for data_type in ['sc']:
            print(f"Calculating metrics for {data_type} data...", flush=True)
            par['prediction'] = prediction_file_name(dataset, data_type)
            rr_store = []
            metric_reg2 = main_reg2(par)
            rr_store.append(metric_reg2)
            _, metric_ws = main_ws_distance(par)
            rr_store.append(metric_ws)
            rr = pd.concat(rr_store, axis=1)
            rr['data_type'] = data_type
            rr_all_store.append(rr)
        rr_all = pd.concat(rr_all_store, axis=0)
        rr_all['dataset'] = dataset
        rr_all.to_csv(f'{results_dir}/metrics_{dataset}.csv', index=False)
        metrics_all.append(rr_all)
    metrics_all = pd.concat(metrics_all, axis=0)
    metrics_all.to_csv(f'{results_dir}/metrics_all.csv', index=False)
        