#!/usr/bin/env python3
import pandas as pd
import numpy as np
import mudata as mu
import anndata as ad
import dictys
import os
import sys
import argparse
import gzip
import subprocess
from pathlib import Path
import warnings
import gdown
import os
import sys
warnings.filterwarnings("ignore")

## VIASH START
par = {
    'rna': 'resources_test/grn_benchmark/inference_data/op_rna.h5ad',
    'atac': 'resources_test/grn_benchmark/inference_data/op_atac.h5ad',
    'temp_dir': 'output/d_test/data',
    'qc': False,
    'prediction': 'output/d_test/grn.h5ad'
}
## VIASH END

try:
    sys.path.append(meta['resources_dir'])
except:
    meta = {
        'resources_dir': './',
        'utils_dir': 'src/utils/'
        }
    sys.path.append(meta['resources_dir'])
    sys.path.append(meta['utils_dir'])
from helper import create_mudata, extract_data, create_p2g_links, prepare_for_model, infer_grn, format_tss
from util import download_annotation



def main(par):
    os.makedirs(par['temp_dir'], exist_ok=True)
    par['annotation_file'] = f"{par['temp_dir']}/gene_annotation.gtf.gz"
    par['tss_file'] = os.path.join(par['temp_dir'], "gene_annotation.tss")
    par['motif_file'] = f"{par['temp_dir']}/motif_file.motif"
    par['mudata'] = os.path.join(par['temp_dir'], "mudata.h5mu")
    # par['p2g'] = os.path.join(par['temp_dir'], "peak2gene.csv")
    par['exp_path'] = os.path.join(par['temp_dir'], "expression.tsv.gz")
    par['pks_path'] = os.path.join(par['temp_dir'], "peaks.bed")
    
    # if True:
    #     # print('Downloading required files...', flush=True)
    #     download_annotation(par)
    #     gdown.download("https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_homer_format_0.0001.motif", par['motif_file'], quiet=False)
    #     print('Formatting TSS file ...', flush=True)
    #     format_tss(par)
    #     print('Creating mudata ...', flush=True)
    #     create_mudata(par['rna'], par['atac'], par['mudata'])

    #     print('extract_data ...', flush=True)
    #     extract_data(par)

    print('create_p2g_links ...', flush=True)    
    create_p2g_links(par)
    # print('prepare_for_model ...', flush=True)
    # rna_output, peaks_output, tfb_path = prepare_for_model(preprocessed_path, p2g_path, par['temp_dir'], par['motif_file'])
    # print('infer_grn ...', flush=True)
    # grn_path = infer_grn(rna_output, peaks_output, tfb_path, par['prediction'])
    # print(f"GRN inference complete. Results saved to {grn_path}")
def create_subsets(par):
    adata = ad.read_h5ad(par['rna'], backed='r')
    # Get all barcodes (cell names) from adata
    barcodes = adata.obs_names.to_list()

    # Write one Subset containing all cells
    subset_dir = f"{par['temp_dir']}/data/subsets/SubsetAll"
    os.makedirs(subset_dir, exist_ok=True)

    with open(f"{subset_dir}/names_rna.txt", "w") as f:
        f.write("\n".join(barcodes))
    
    import shutil
    shutil.copy(f"{subset_dir}/names_rna.txt", f"{subset_dir}/names_atac.txt")

    # Also write subsets.txt (only one entry: SubsetAll)
    with open(f"{par['temp_dir']}/data/subsets.txt", "w") as f:
        f.write("SubsetAll\n")
if __name__ == "__main__":
    main(par)
    # python -m dictys chromatin tssdist output/d_test/data/expression.tsv.gz output/d_test/data/peaks.bed output/d_test/data/gene_annotation.tss output/d_test/data/tssdist.tsv.gz
    