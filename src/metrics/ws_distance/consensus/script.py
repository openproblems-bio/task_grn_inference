import pandas as pd 
import anndata as ad 
from tqdm import tqdm
import numpy as np
import os


import argparse

naming_convention = lambda dataset, method: f'{dataset}.{method}.{method}.prediction.h5ad'

arg = argparse.ArgumentParser(description='Compute consensus number of putative regulators for each gene')
arg.add_argument('--dataset', type=str, help='Dataset to use for the analysis')
arg.add_argument('--models_dir', type=str, help='Directory containing the GRN models')
arg.add_argument('--evaluation_data_sc', type=str, help='Path to the evaluation data')
arg.add_argument('--ws_consensus', type=str, help='Path to save the consensus regulators')
arg.add_argument('--tf_all', type=str, help='Path to the file containing all transcription factors')
arg.add_argument('--models', nargs='+', help='List of models to use for the analysis')
args = arg.parse_args()

par = args.__dict__

def main(par):
    adata = ad.read_h5ad(par['evaluation_data_sc'])
    tf_all = np.loadtxt(par['tf_all'], dtype='str')
    available_tfs = np.intersect1d(adata.obs['perturbation'].unique(), tf_all)

    # - read all models
    grn_store = []
    for model in par['models']:
        prediction_file = f"{par['models_dir']}/{naming_convention(par['dataset'], model)}"
        if not os.path.exists(prediction_file):
            print('Warnings: ', prediction_file, ' doesnt exists')
            continue
        else:
            grn = ad.read_h5ad(prediction_file)
            grn = pd.DataFrame(grn.uns['prediction'])

        grn['model'] = model
        grn_store.append(grn)
    grn_all = pd.concat(grn_store).reset_index(drop=True)

    assert len(grn_all) > 0, 'No GRN predictions found in the models'
    # - subset to available TFs
    grn_all = grn_all[grn_all['source'].isin(available_tfs)]

    # - determine consensus
    edges_count = grn_all.groupby(['source', 'model']).size().reset_index(name='n_edges').pivot(index='source',columns='model').fillna(0)

    consensus = []
    for tf, row in edges_count.iterrows():
        row_nozero = row[row!=0] 
        consensus.append({'source':tf, 'theta':'ws-theta-0.0', 'value':int(np.quantile(row_nozero, 0))})
        consensus.append({'source':tf, 'theta':'ws-theta-0.5', 'value':int(np.quantile(row_nozero, 0.5))})
        consensus.append({'source':tf, 'theta':'ws-theta-1.0', 'value':int(np.quantile(row_nozero, 1))})
    consensus = pd.DataFrame(consensus)
    # - save 
    print('saving the consensus to ', par['ws_consensus'])
    consensus.to_csv(par['ws_consensus'])

    return consensus

if __name__ == '__main__':
    main(par)
