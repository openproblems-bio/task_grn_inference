import pereggrn_perturbations
import pereggrn_networks

import os 
import anndata as ad 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import scanpy as sc

pereggrn_dir = '/home/jnourisa/projs/ongoing/pereggrn'

pereggrn_perturbations.set_data_path(f"{pereggrn_dir}/perturbation_data/perturbations")

pereggrn_perturbations.load_perturbation_metadata()

meta = {
    "resources_dir": 'src/utils/',
}
sys.path.append(meta["resources_dir"])
from util import process_links

def get_dataset(par):
    for file_name in par['datasets']:

        adata = pereggrn_perturbations.load_perturbation(file_name) 
        # pereggrn_perturbations.check_perturbation_dataset(ad = adata)
        adata.write(f"{par['raw_datasets_dir']}/{file_name}.h5ad")

def get_networks(par):
    os.makedirs(par['nets_dir'], exist_ok=True)

    names = []
    for model in par['nets']:
        net = pd.read_parquet(f"{pereggrn_dir}/network_collection/networks/{model}")
        net.columns = ['source','target','weight']
        method = model.split('/')[0].split('_')[0].capitalize()
        tissue = model.split('/')[-1].split('.')[0].replace('_', ' ').capitalize()
        name = method+':'+tissue

        net = process_links(net, par)

        net.to_csv(f"{par['nets_dir']}/{name}.csv")



def main(par):
    print('Getting data ...')
    get_dataset(par)
    print('Getting networks ...')
    get_networks(par)
        
if __name__ == '__main__':
    par = {
        'datasets': ['norman', 'adamson', 'nakatake'],
        'nets': [
                            'ANANSE_tissue/networks/lung.parquet',
                            'ANANSE_tissue/networks/stomach.parquet', 
                            'ANANSE_tissue/networks/heart.parquet',
                            'ANANSE_tissue/networks/bone_marrow.parquet',
                            
                            'gtex_rna/networks/Whole_Blood.parquet',
                            'gtex_rna/networks/Brain_Amygdala.parquet', 
                            'gtex_rna/networks/Breast_Mammary_Tissue.parquet', 
                            'gtex_rna/networks/Lung.parquet',
                            'gtex_rna/networks/Stomach.parquet',

                            'cellnet_human_Hg1332/networks/bcell.parquet',
                            'cellnet_human_Hg1332/networks/tcell.parquet',
                            'cellnet_human_Hg1332/networks/skin.parquet',
                            'cellnet_human_Hg1332/networks/neuron.parquet',
                            'cellnet_human_Hg1332/networks/heart.parquet',
                            ],
        'raw_datasets_dir': 'resources/datasets_raw/',
        'nets_dir': f'resources/grn_models/global/',
        'max_n_links': 50_000,
    }
    
    

    
    main(par)