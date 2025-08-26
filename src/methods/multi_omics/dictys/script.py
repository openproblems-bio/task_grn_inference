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
    'rna': 'resources/grn_benchmark/inference_data/op_rna.h5ad',
    'atac': 'resources/grn_benchmark/inference_data/op_atac.h5ad',
    'distance': 500000,
    'temp_dir': 'output/d_test/',
    'qc': False,
    'prediction': 'output/d_test/grn.h5ad'
}
## VIASH END

try:
    sys.path.append(meta['resources_dir'])
except:
    meta = {
        'resources_dir': './'
        }
    sys.path.append(meta['resources_dir'])
from helper import create_mudata, preprocess_data, extract_data, create_p2g_links, prepare_for_model, infer_grn
# from src.methods.multi_omics.dictys.helper import create_mudata, preprocess_data, extract_data, create_p2g_links, prepare_for_model, infer_grn
def main(par):
    
    output_dir = par['temp_dir']
    os.makedirs(output_dir, exist_ok=True)
    par['gene_annotation'] = f"{output_dir}/gene_annotation.gtf.gz"
    par['tss_file'] = os.path.join(output_dir, "gene_annotation.tss")
    par['motif_file'] = f"{output_dir}/motif_file.motif"
    par['mudata'] = os.path.join(output_dir, "mudata.h5mu")
    par['p2g'] = os.path.join(output_dir, "peak2gene.csv")
    # par['preprocessed_mudata'] = os.path.join(output_dir, "preprocessed.h5mu")
    par['exp_path'] = os.path.join(output_dir, "expression.tsv.gz")
    par['pks_path'] = os.path.join(output_dir, "peaks.bed")
    

    cmd = f"""zcat {par['gene_annotation']} \
    | awk '$3 == "transcript" {{
        split($9,a,";");
        for (i in a) {{ if (a[i] ~ /gene_id/) {{ gsub(/"/,"",a[i]); gsub(/gene_id /,"",a[i]); gene=a[i] }} }}
        if ($7 == "+") {{tss=$4}} else {{tss=$5}}
        print gene "\\t" $1 "\\t" tss "\\t" $7
    }}' > {par['tss_file']}"""

    df_tss = pd.read_csv(
        par['tss_file'], 
        sep='\t', 
        header=None, 
        usecols=[0, 1, 2, 3],  # columns: gene, chrom, TSS, strand
        names=['gene', 'chrom', 'tss', 'strand']
    )

    # Keep only rows where chromosome is numeric
    df_tss = df_tss[df_tss['chrom'].str.isnumeric()]

    # Prepend "chr" to the chromosome column
    df_tss['chrom'] = "chr" + df_tss['chrom']

    df_tss.reset_index(drop=True, inplace=True)

    df_tss.to_csv(par['tss_file'], sep='\t', index=False, header=False)


    # run it
    os.system(cmd)
    
    if True:
        print('Downloading required files...', flush=True)
        # gdown.download("http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz", par['gene_annotation'], quiet=False)
        # gdown.download("https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_homer_format_0.0001.motif", par['motif_file'], quiet=False)

        

        # Step 1: Create mudata from RNA and ATAC
        # print('Creating mudata ...', flush=True)
        # create_mudata(par['rna'], par['atac'], par['mudata'])

        # Step 2: Preprocess data
        # print('preprocess_data ...', flush=True)
        # preprocess_data(mudata_path, output_dir, par['qc'])

        # Step 3: Extract data
        print('extract_data ...', flush=True)
        extract_data(par)
    # Step 4: Create peak-to-gene linkages
    # --- here
    print('create_p2g_links ...', flush=True)    
    p2g_path = create_p2g_links(par)

    # Step 5: Prepare data for modeling
    print('prepare_for_model ...', flush=True)
    rna_output, peaks_output, tfb_path = prepare_for_model(preprocessed_path, p2g_path, output_dir, par['motif_file'])

    # Step 6: Infer GRN
    print('infer_grn ...', flush=True)
    grn_path = infer_grn(rna_output, peaks_output, tfb_path, par['prediction'])

    print(f"GRN inference complete. Results saved to {grn_path}")
        
if __name__ == "__main__":
    main(par)