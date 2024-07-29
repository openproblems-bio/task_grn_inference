import sys
import anndata as ad
import scanpy as sc
import os 
import pandas as pd
import numpy as np
from scipy.io import mmwrite

import argparse


par = {
  "multiomics_rna": "resources_test/grn-benchmark/multiomics_rna.h5ad",
  "multiomics_atac": "resources_test/grn-benchmark/multiomics_atac.h5ad",
  "rna_matrix": 'output/scRNA/X_matrix.mtx',
  "atac_matrix": 'output/scATAC/X_matrix.mtx',
  "rna_gene_annot": 'output/scRNA/annotation_gene.csv',
  "rna_cell_annot": 'output/scRNA/annotation_cell.csv',
  "atac_peak_annot": 'output/scATAC/annotation_gene.csv',
  "atac_cell_annot": 'output/scATAC/annotation_cell.csv'

}
print(par)
def format_data(par):
    os.makedirs("output/scATAC/", exist_ok=True)
    os.makedirs("output/scRNA/", exist_ok=True)
    
    print('Reading input files', flush=True)
    rna = ad.read_h5ad(par['multiomics_rna'])
    atac = ad.read_h5ad(par['multiomics_atac'])

    # save sparse matrix
    print(atac)
    mmwrite(f"{par['atac_matrix']}", atac.X)
    # save annotation
    annotation_peak = atac.var.reset_index().location.str.split(':', expand=True)
    annotation_peak.columns = ['seqname', 'ranges']
    annotation_peak['strand'] = '+' 
    annotation_peak.to_csv(f"{par['atac_peak_annot']}")

    annotation_cells = atac.obs.reset_index()
    annotation_cells.to_csv(f"{par['atac_cell_annot']}")


    # save sparse matrix
    print(rna)
    mmwrite(f"{par['rna_matrix']}", rna.X)
    # save annotation
    annotation_gene = rna.var.reset_index()
    annotation_gene.to_csv(f"{par['rna_gene_annot']}")


    annotation_cells = rna.obs.reset_index()[['obs_id','cell_type']]
    annotation_cells.to_csv(f"{par['rna_cell_annot']}")

    print('Format data completed', flush=True)

format_data(par)


