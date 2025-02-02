import lamindb as ln
ln.setup.init(storage='./output')

from scprint import scPrint
from scdataloader import Preprocessor, utils
from scprint.tasks import GNInfer, Embedder, Denoiser, withknn
from bengrn import BenGRN, get_sroy_gt, compute_genie3
from bengrn.base import train_classifier
from grnndata import utils as grnutils
from grnndata import read_h5ad
from anndata.utils import make_index_unique
from anndata import concat
import scanpy as sc
from matplotlib import pyplot as plt
import numpy as np
import anndata as ad
import pandas as pd
from scdataloader.utils import populate_my_ontology

import torch
torch.set_float32_matmul_precision('medium')
import sys 
sys.path.append("./")
from src.helper_infer_grns import efficient_melting 
## VIASH START
par = {
    'rna': '../task_grn_inference/resources/inference_datasets/op_rna.h5ad',
    'tf_all': '../task_grn_inference/resources/prior/tf_all.csv',
    'prediction': 'output/grn.h5ad',
    'filtration': 'top-k',
    'max_num_links': 50_000,
    'num_genes': 5000,
    'max_cells': 1000
}
## VIASH END
def run_scprint_sub(adata, model, par):
    grn_inferer = GNInfer(
                        how="most var within",
                        preprocess="softmax",
                        head_agg='mean',
                        forward_mode="none",
                        filtration=par['filtration'],
                        num_genes=par['num_genes'],
                        num_workers=8, 
                        max_cells=par['max_cells'],
                        doplot=False,
                        batch_size=16,
                        k=par['max_num_links'],
                        )

    grn = grn_inferer(model, adata)

    net = grn.varp['GRN'].todense().A
    print(net.shape)
    # - melt to have the format of the benchmark
    gene_names = grn.var['gene_name'].values
    net = efficient_melting(net, gene_names)
    # - subset to TFs
    tfs = np.loadtxt(par["tf_all"], dtype=str)
    tf_names = [gene_name for gene_name in gene_names if (gene_name in tfs)]
    net = net[net['source'].isin(tf_names)]
    return net
def run_scprint(par):
    adata = ad.read_h5ad(par['rna'])
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

    preprocessor = Preprocessor(do_postp=False, is_symbol=True)
    adata = preprocessor(adata)

    model = scPrint.load_from_checkpoint('output/medium.ckpt', 
    precpt_gene_emb = None)

    for i, cell_type in enumerate(adata.obs["cell_type"].unique()):
        print(cell_type)
        adata_cell_type = adata[adata.obs["cell_type"] == cell_type].copy()
        net = run_scprint_sub(adata_cell_type, model, par)
        net['cell_type'] = cell_type
        if i == 0:
            net_all = net
        else:
            net_all = pd.concat([net_all, net], axis=0)

    return net_all 

    
if __name__ == '__main__':
    populate_my_ontology( 
        organisms= ["NCBITaxon:9606"]
    )
    
    net = run_scprint(par)

    net['weight'] = net['weight'].astype(str)
    output = ad.AnnData(X=None, uns={"method_id": par['method_id'], "dataset_id": par['dataset_id'], "prediction": net[["source", "target", "weight"]]})
    output.write(par['prediction'])