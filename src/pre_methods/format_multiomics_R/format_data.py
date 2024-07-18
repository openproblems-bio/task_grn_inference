import sys
import anndata as ad
import scanpy as sc
import os 
import pandas as pd
import numpy as np
from scipy.io import mmwrite

par = {
  "multiomics_rna": "resources/grn-benchmark/multiomics_rna.h5ad",
  "multiomics_atac": "resources/grn-benchmark/multiomics_atac.h5ad",
  "temp_dir": 'output/figr'

}

def format_data(par):
    os.makedirs(f"{par['temp_dir']}/scATAC/", exist_ok=True)
    os.makedirs(f"{par['temp_dir']}/scRNA/", exist_ok=True)
    
    os.makedirs(par['temp_dir'], exist_ok=True)
    print('Reading input files', flush=True)
    rna = ad.read_h5ad(par['multiomics_rna'])
    atac = ad.read_h5ad(par['multiomics_atac'])

    # save sparse matrix
    mmwrite(f"{par['temp_dir']}/scATAC/X_matrix.mtx", atac.X)
    # save annotation
    annotation_peak = atac.var.reset_index().location.str.split(':', expand=True)
    annotation_peak.columns = ['seqname', 'ranges']
    annotation_peak['strand'] = '+' 
    annotation_peak.to_csv(f"{par['temp_dir']}/scATAC/annotation_peak.csv")

    annotation_cells = atac.obs.reset_index()
    annotation_cells.to_csv(f"{par['temp_dir']}/scATAC/annotation_cells.csv")


    # save sparse matrix
    mmwrite(f"{par['temp_dir']}/scRNA/X_matrix.mtx", rna.X)
    # save annotation
    annotation_gene = rna.var.reset_index()
    annotation_gene.to_csv(f"{par['temp_dir']}/scRNA/annotation_gene.csv")


    annotation_cells = rna.obs.reset_index()[['obs_id','cell_type']]
    annotation_cells.to_csv(f"{par['temp_dir']}/scRNA/annotation_cells.csv")

    print('Format data completed', flush=True)

format_data(par)