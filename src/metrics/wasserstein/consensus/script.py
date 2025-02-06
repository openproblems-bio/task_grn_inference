import pandas as pd 
import anndata as ad 
from tqdm import tqdm
import numpy as np
import os

# - determine consensus 
par = {
    'evaluation_data_sc': f'resources/datasets_raw/norman_sc_counts.h5ad',
    'ws_consensus': 'resources/grn_benchmark/prior/consensus_ws_distance_norman.csv',
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
    'models_dir': 'resources/grn_models/norman',
    'models': ['pearson_corr', 'grnboost2','portia', 'ppcor','scenic']
}
def main(par):

    adata = ad.read_h5ad(par['evaluation_data_sc'])
    tf_all = np.loadtxt(par['tf_all'], dtype='str')
    available_tfs = np.intersect1d(adata.obs['perturbation'].unique(), tf_all)

    # - read all models
    grn_store = []
    for model in par['models']:
        prediction_file = f"{par['models_dir']}/{model}.csv"
        if not os.path.exists(prediction_file):
            print(prediction_file, ' doesnt exists')
            continue
        else:
            grn = pd.read_csv(prediction_file, index_col=0)

        grn['model'] = model
        grn_store.append(grn)
    grn_all = pd.concat(grn_store).reset_index(drop=True)
    # print(grn_all)

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
    consensus.to_csv(par['ws_consensus'])

    return consensus

if __name__ == '__main__':
    main(par)