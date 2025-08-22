import os
import pandas as pd
import numpy as np
import sys
import anndata as ad
import scanpy as sc
from tqdm import tqdm




meta = {
    "resources_dir": './',
}
sys.path.append(meta["resources_dir"])

from src.control_methods.pearson_corr.script import main as main_inference
from src.metrics.regression_2.helper import main as main_reg2
from src.metrics.ws_distance.helper import main as main_ws_distance
from src.metrics.ws_distance.consensus.helper import main as main_consensus_ws_distance
from src.process_data.helper_data import sum_by

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
def prediction_file_name(dataset, granularity):
    return f'{results_dir}/{dataset}.prediction_{granularity}.h5ad'

def main_pseudobulk(par):

    # - read inputs and cluster with differen resolutions
    rna = ad.read_h5ad(par['rna'])
    granularity = par['granularity']

    if granularity == -1:
        pass
    else:
        sc.pp.pca(rna, layer=par['layer'])
        sc.pp.neighbors(rna)
        sc.tl.leiden(rna, resolution=granularity, key_added=f'leiden_{granularity}')
        rna_bulk = sum_by(rna, f'leiden_{granularity}', unique_mapping=True)
        rna_bulk.layers[par['layer']] = rna_bulk.X
        for key in rna.uns.keys():
            rna_bulk.uns[key] = rna.uns[key]
        rna = rna_bulk

    rna.write(par['rna_pseudobulked'])


if __name__ == '__main__':
    results_dir = 'resources/results/experiment/granular_pseudobulk'
    os.makedirs(results_dir, exist_ok=True)
    dataset = 'op'

    par = def_par(dataset)
    degrees = [-1.0, 1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0]
    # - pseudobulk
    if False:
        print('Pseudobulking data...', flush=True)
        par['rna'] = f'resources/grn_benchmark/inference_data/{dataset}_rna.h5ad'
        print(ad.read_h5ad(par['rna']))
        
        for granuality in tqdm(degrees, desc='Pseudobulking'):
            par['granularity'] = granuality
            par['rna_pseudobulked'] = f'{results_dir}/{dataset}_granularity_{granuality}.h5ad'
            main_pseudobulk(par)
    # - infer grns
    if False:
        print('Inferring GRNs for pseudobulked data...', flush=True)
        for granuality in tqdm(degrees, desc='Inferring GRNs'):
            par['rna'] = f'{results_dir}/{dataset}_granularity_{granuality}.h5ad'
            net = main_inference(par)

            par['prediction'] = prediction_file_name(dataset, granuality)
            net.write_h5ad(par['prediction'])
    
    # - consensus 
    par['regulators_consensus'] = f'{results_dir}/regulators_consensus_{dataset}.json'
    if False:
        print('Calculating consensus for pseudobulked data...', flush=True)
        def naming_convention(dataset, model):
            return f'{dataset}.{model}.h5ad'
        from src.metrics.regression_2.consensus.helper import main as main_consensus_reg2
        par['naming_convention'] = naming_convention
        par['dataset'] = dataset
        par['models_dir'] = results_dir
        par['models'] = [f'prediction_{i}' for i in degrees]
        
        _ = main_consensus_reg2(par)
    
        # par['ws_consensus'] = f'{results_dir}/ws_consensus_{dataset}.json'
        # main_consensus_ws_distance(par)

    # - grn evaluation
    if True:
        print('Calculating metrics for pseudobulked data...', flush=True)
        rr_all_store = []
        for granuality in tqdm(degrees, desc='Calculating metrics'):
            print(f"Calculating metrics for {granuality} ...")
            par['prediction'] = prediction_file_name(dataset, granuality)
            metric_reg2 = main_reg2(par)
            # _, metric_ws = main_ws_distance(par)

            rr_store = []
            rr_store.append(metric_reg2)
            # rr_store.append(metric_ws)
            rr = pd.concat(rr_store, axis=1)
            rr['granularity'] = granuality
            rr_all_store.append(rr)
        rr_all = pd.concat(rr_all_store, axis=0)
        rr_all['dataset'] = dataset
        rr_all.to_csv(f'{results_dir}/metrics_{dataset}.csv', index=False)
