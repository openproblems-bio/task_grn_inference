import sys
import anndata as ad
import scanpy as sc
import os 
import pandas as pd
import numpy as np
from scipy.io import mmwrite

import argparse

# Initialize the argument parser
parser = argparse.ArgumentParser(description="Update dictionary with command-line arguments.")

# Define the expected arguments
parser.add_argument("--multiomics_rna", required=True, help="Path to the multiomics RNA file")
parser.add_argument("--multiomics_atac", required=True, help="Path to the multiomics ATAC file")
parser.add_argument("--rna_matrix", required=True, help="Path to the RNA matrix file")
parser.add_argument("--atac_matrix", required=True, help="Path to the ATAC matrix file")
parser.add_argument("--rna_gene_annot", required=True, help="Path to the RNA gene annotation file")
parser.add_argument("--rna_cell_annot", required=True, help="Path to the RNA cell annotation file")
parser.add_argument("--atac_peak_annot", required=True, help="Path to the ATAC peak annotation file")
parser.add_argument("--atac_cell_annot", required=True, help="Path to the ATAC cell annotation file")

# Parse the arguments
args = parser.parse_args()

# Create the dictionary with the parsed arguments
par = {
    "multiomics_rna": args.multiomics_rna,
    "multiomics_atac": args.multiomics_atac,
    "rna_matrix": args.rna_matrix,
    "atac_matrix": args.atac_matrix,
    "rna_gene_annot": args.rna_gene_annot,
    "rna_cell_annot": args.rna_cell_annot,
    "atac_peak_annot": args.atac_peak_annot,
    "atac_cell_annot": args.atac_cell_annot
}

# par = {
#   "multiomics_rna": "resources/grn-benchmark/multiomics_rna.h5ad",
#   "multiomics_atac": "resources/grn-benchmark/multiomics_atac.h5ad",
#   "rna_matrix": 'output/scRNA/X_matrix.mtx',
#   "atac_matrix": 'output/scATAC/X_matrix.mtx',
#   "rna_gene_annot": 'output/scRNA/annotation_gene.csv',
#   "rna_cell_annot": 'output/scRNA/annotation_cell.csv',
#   "atac_peak_annot": 'output/scATAC/annotation_gene.csv',
#   "atac_cell_annot": 'output/scATAC/annotation_cell.csv'

# }
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


