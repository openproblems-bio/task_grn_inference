import os
import sys 
from anndata.utils import make_index_unique
from anndata import concat
import scanpy as sc
from matplotlib import pyplot as plt
import numpy as np
import anndata as ad
import pandas as pd
import requests
import torch
import argparse


import lamindb_setup as ln_setup
ln_setup.init(storage="./testdb", name="test", modules="bionty")

from scdataloader.utils import populate_my_ontology
from scprint import scPrint
from scdataloader import Preprocessor, utils
from scprint.tasks import GNInfer, Embedder, Denoiser, withknn
from bengrn import BenGRN, get_sroy_gt, compute_genie3
from bengrn.base import train_classifier
from grnndata import utils as grnutils
from grnndata import read_h5ad

torch.set_float32_matmul_precision('medium')


## VIASH START
par = {
    'rna': 'resources/grn_benchmark/inference_data/op_rna.h5ad',
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
    'prediction': 'output/grn.h5ad',
    'filtration': 'none', #'top-k',
    'max_n_links': 50_000,
    'num_genes': 5000,
    'max_cells': 1000,
    'num_workers': 20,
    'model': 'medium',
    'download_checkpoint': False,
    'populate_ontology': False,
    'temp_dir': 'output/scprint/',
    'method_id': 'scprint',
    'dataset_id': 'op'
}
## VIASH END

## LOCAL START
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=int, help='Number of cores')
parser.add_argument('--max_n_links', type=int, help='Number of top links to retain')
parser.add_argument('--dataset_id', type=str, help='Dataset id')
parser.add_argument('--normalize', action='store_true')
args = parser.parse_args()

par_local = vars(args)

for key, value in par_local.items():
    if value is not None:
        par[key] = value

## LOCAL END

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
    'resources_dir': 'src/utils'
    }
    sys.path.append(meta["resources_dir"])
from util import process_links, efficient_melting


def main_sub(adata, model, par):
    grn_inferer = GNInfer(
                        how="most var within",
                        preprocess="softmax",
                        head_agg='mean',
                        forward_mode="none",
                        filtration=par['filtration'],
                        num_genes=par['num_genes'],
                        num_workers=par['num_workers'], 
                        max_cells=par['max_cells'],
                        doplot=False,
                        batch_size=16,
                        precision="32-mixed",
                        )

    grn = grn_inferer(model, adata)

    net = grn.varp['GRN']
    if hasattr(net, "todense"):  # Check if it's a sparse matrix
        net = net.todense().A.T
    else:  # Already a NumPy array
        net = net.T
    
    # - melt to have the format of the benchmark
    gene_names = grn.var['gene_name'].values
    net = efficient_melting(net, gene_names)
    net = net[net['weight'] != 0]
    assert ~net[['source', 'target', 'weight']].duplicated().any()
    

    # - subset to TFs
    tfs = np.loadtxt(par["tf_all"], dtype=str)
    tf_names = [gene_name for gene_name in gene_names if (gene_name in tfs)]
    net = net[net['source'].isin(tf_names)]

    net = net.sort_values(by='weight', ascending=False, key=abs)[:par['max_n_links']]
    
    return net
def main(par):
    adata = ad.read_h5ad(par['rna'])
    if hasattr(adata, "layers"):
        if 'counts' in adata.layers:
            adata.X = adata.layers['counts']
    # if hasattr(adata, "raw"):
    #     del adata.raw
    adata.var["gene_name"] = adata.var.index

    adata.obs['organism_ontology_term_id'] = 'NCBITaxon:9606'
    # adata.obs['self_reported_ethnicity_ontology_term_id'] = "HANCESTRO:0005"
    # adata.obs['sex_ontology_term_id'] = "PATO:0000384"
    # adata.obs['disease_ontology_term_id'] = "MONDO:0000001"
    # adata.obs['development_stage_ontology_term_id'] = "HsapDv:0000087"
    # adata.obs['tissue_ontology_term_id'] = "UBERON:0000178"
    # adata.obs['assay_ontology_term_id'] = "unknown"


    # cell_type_to_ontology = {
    #     "T cells": "CL:0000084",
    #     "B cells": "CL:0000236",
    #     "Myeloid cells": "CL:0000763",
    #     "NK cells": "CL:0000623",
    # }

    # adata.obs["cell_type_ontology_term_id"] = adata.obs["cell_type"].apply(lambda name: cell_type_to_ontology.get(name, name))
    min_valid_genes = min(adata.X.shape[1]-1, 2000)
    print(f'Min valid genes: {min_valid_genes} ')
    preprocessor = Preprocessor(do_postp=False, is_symbol=True, skip_validate=True, 
                                force_preprocess=False, use_raw=False, min_valid_genes_id=min_valid_genes, min_dataset_size=1024)
    adata = preprocessor(adata)

    model = scPrint.load_from_checkpoint(par['checkpoint'], 
    precpt_gene_emb = None)

    if 'cell_type' not in adata.obs:
        adata.obs['cell_type'] = 'dummy_cell_type'
    for i, cell_type in enumerate(adata.obs["cell_type"].unique()):
        print(cell_type)
        adata_cell_type = adata[adata.obs["cell_type"] == cell_type].copy()
        net = main_sub(adata_cell_type, model, par)
        net['cell_type'] = cell_type
        if i == 0:
            net_all = net
        else:
            net_all = pd.concat([net_all, net], axis=0)

    return net_all 

    
if __name__ == '__main__':
    dataset_id = ad.read_h5ad(par['rna'], backed='r').uns['dataset_id']
    if par['populate_ontology']: 
        populate_my_ontology( 
            organisms= ["NCBITaxon:9606"]
        )
    par['checkpoint'] = par['temp_dir'] + '/scprint.ckpt'
    if par['download_checkpoint']:

        os.makedirs(par['temp_dir'], exist_ok=True)
        print(f"Downloading checkpoint")
        checkpoint_link = f"https://huggingface.co/jkobject/scPRINT/resolve/main/{par['model']}.ckpt" #TODO: experiment with this
        
        response = requests.get(checkpoint_link, stream=True)
        if response.status_code == 200:
            with open(par['checkpoint'], 'wb') as f:
                for chunk in response.iter_content(chunk_size=1024):
                    f.write(chunk)
            print(f"Downloaded checkpoint to {par['checkpoint']}")
        else:
            print(f"Error: Failed to download checkpoint (status code {response.status_code})")

    net = main(par)

    print(f"Writing results to {par['prediction']}")
    net = net.astype(str)
    output = ad.AnnData(X=None, uns={"method_id": "scprint", "dataset_id": dataset_id, "prediction": net[["source", "target", "weight", "cell_type"]]})
    output.write(par['prediction'])