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
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=UserWarning, module="torch")
def format_tss(par):
    import gzip
    seen = set()  # to track unique gene names
    with gzip.open(par['annotation_file'], 'rt') as f_in, open(par['tss_file'], 'w') as f_out:
        for line in f_in:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if parts[2] != "gene":
                continue

            attr_str = parts[8]
            gene_name = ""
            for attr in attr_str.split(";"):
                attr = attr.strip()
                if attr.startswith('gene_name'):
                    gene_name = attr.split('"')[1]
                    break
            if not gene_name:  # fallback to gene_id
                for attr in attr_str.split(";"):
                    attr = attr.strip()
                    if attr.startswith('gene_id'):
                        gene_name = attr.split('"')[1]
                        break

            # skip duplicates
            if gene_name in seen:
                continue
            seen.add(gene_name)

            chrom = parts[0]
            # Add "chr" if chromosome is numeric, X, or Y
            if chrom.isdigit() or chrom in ["X", "Y"]:
                chrom = "chr" + chrom

            start = int(parts[3]) - 1
            end = parts[4]
            strand = parts[6]

            f_out.write(f"{chrom}\t{start}\t{end}\t{gene_name}\t.\t{strand}\n")

    print(f"TSS file written to {par['tss_file']} with {len(seen)} unique genes")
def create_mudata(rna_path, atac_path, output_path):
    """Combine RNA and ATAC data into mudata object"""
    print("Creating mudata object from RNA and ATAC data...")
    
    # Read input files
    rna = ad.read_h5ad(rna_path)
    atac = ad.read_h5ad(atac_path)
    print(f"RNA shape: {rna.shape}, ATAC shape: {atac.shape}", flush=True)
    
    # Ensure counts layer exists
    if 'counts' not in rna.layers:
        print("Creating counts layer in RNA data...")
        rna.layers['counts'] = rna.X.copy()
    
    if 'counts' not in atac.layers:
        print("Creating counts layer in ATAC data...")
        atac.layers['counts'] = atac.X.copy()
    
    # Create mudata object
    mdata = mu.MuData({'rna': rna, 'atac': atac})
    
    # Write to file
    mdata.write(output_path)
    
    
def extract_data(par, use_p2g=False, p2g_path=None):
    """Extract RNA and ATAC data for processing"""
    print("Extracting RNA and ATAC data...")

    # Read TSS annotation
    tss_genes = set(pd.read_csv(par['tss_file'], header=None, index_col=None, sep='\t')[3])

    # Write the RNA matrix
    data = mu.read(par['mudata'])
    rna = data['rna']
    rna.var_names_make_unique()
    rna_X = pd.DataFrame(
        np.array(rna.layers['counts'].todense()).T,
        columns=data['rna'].obs.index,
        index=data['rna'].var.index
    )

    # Filter RNA to keep only genes present in TSS annotation
    
    rna_X = rna_X.loc[rna_X.index.isin(tss_genes)]
    print(f"RNA data shape after filtering with TSS genes: {rna_X.shape}", flush=True)
    rna_X.to_csv(par['exp_path'], sep="\t", compression="gzip")

    # ATAC peaks
    if use_p2g and p2g_path:
        all_atac_peak = np.unique(pd.read_csv(p2g_path)['cre'])
    else:
        all_atac_peak = np.unique([n.replace(':', '-') for n in data['atac'].var.index])

    all_atac_peak = pd.DataFrame([n.split('-') for n in all_atac_peak])
    all_atac_peak.columns = ['chr', 'srt', 'end']
    all_atac_peak['srt'] = all_atac_peak['srt'].astype(int)
    all_atac_peak['end'] = all_atac_peak['end'].astype(int)
    all_atac_peak = all_atac_peak[(all_atac_peak.end - all_atac_peak.srt) >= 100]
    all_atac_peak = all_atac_peak.sort_values(by=['chr', 'srt', 'end'])
    print(f"Peaks: {len(all_atac_peak)}", flush=True)
    all_atac_peak.to_csv(par['pks_path'], sep='\t', header=False, index=False)

    # Store clusters if celltype annotation exists
    if 'celltype' in data.obs.columns:
        clus = sorted(data.obs['celltype'].unique())
        for c in clus:
            ctype_ids = data[data.obs['celltype'] == c].obs.index
            c = c.replace(' ', '_')
            with open(os.path.join(par['temp_dir'], f'barcodes_{c}.txt'), "w") as f:
                for i in ctype_ids:
                    f.write(f"{i}\n")
    
def create_p2g_links(par):
    """Create peak-to-gene linkages based on genomic proximity"""
    print("Creating peak-to-gene linkages...")
    
    atac_filename = os.path.join(par['temp_dir'], "atac_peak.tsv.gz")
    dist_filename = os.path.join(par['temp_dir'], "tssdist.tsv.gz")
    
    # Get ATAC peaks from BED file
    atac_peaks = pd.read_csv(par['pks_path'], sep='\t', header=None)
    atac_peaks.columns = ['chr', 'start', 'end']
    atac_peak_names = [f"{row['chr']}:{row['start']}:{row['end']}" for _, row in atac_peaks.iterrows()]
    
    # Create empty ATAC matrix for dictys
    atac_X = pd.DataFrame(np.zeros((len(atac_peak_names), 1)), index=atac_peak_names, columns=['placeholder'])
    atac_X.to_csv(atac_filename, sep="\t", compression="gzip")
    
    # Run dictys to identify peaks that are within specified distance of annotated TSS
    ret = os.system(
        f'PYTHONWARNINGS="ignore" python3 -m dictys chromatin tssdist '
        f"{par['exp_path']} {atac_filename} {par['tss_file']} {dist_filename} "
    )

    if ret != 0:
        print(f"Command failed with exit code {ret}", file=sys.stderr)
        sys.exit(ret)  # stops the script immediately
    # Convert distance to score for p2g
    df = pd.read_csv(dist_filename, sep='\t').rename(columns={'region': 'cre', 'target': 'gene', 'dist': 'score'})
    df['score'] = -np.abs(df['score'])
    df['cre'] = df['cre'].str.replace(':', '-')
    df = df.sort_values('score', ascending=False).reset_index(drop=True).reset_index(names='rank')
    df['score'] = (1 - (df['rank'] / df['rank'].max()))
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(par['p2g']), exist_ok=True)
    df[['cre', 'gene', 'score']].to_csv(par['p2g'], index=False)
    
def prepare_for_model(preprocessed_path, output_dir, tf_database):
    """Prepare data for modeling"""
    print("Preparing data for GRN inference...")
    
    # Output paths
    rna_output = os.path.join(output_dir, "rna_processed.tsv.gz")
    peaks_output = os.path.join(output_dir, "peaks_processed.tsv.gz")
    tfb_output = os.path.join(output_dir, "tf_binding.tsv")
    
    # Read and process gex
    rna_path = preprocessed_path
    rna = mu.read(rna_path)
    rna.mod['rna'].X = rna.mod['rna'].layers['counts']
    rna_df = rna.mod['rna'].to_df().T
    rna_df.to_csv(rna_output, sep='\t', compression='gzip')
    
    # Process peaks from p2g
    use_peaks = True
    if use_peaks:
        peaks = pd.read_csv(par['p2g'])['cre'].unique()
    else:
        atac_path = preprocessed_path
        peaks = mu.read(atac_path).mod['atac'].var_names
        
    peaks = np.array([p.replace('-', ':') for p in peaks])
    peaks = pd.DataFrame(np.zeros((peaks.size, 1)), index=peaks, columns=['placeholder'])
    peaks.to_csv(peaks_output, sep='\t', compression='gzip')
    
    # Process TF binding information
    tfb = pd.read_csv(tf_database)
    tfb['cre'] = tfb['cre'].str.replace('-', ':')
    tfb = tfb[tfb['tf'].isin(rna_df.index) & tfb['cre'].isin(peaks.index)]
    output_tfb = tfb.rename(columns={'cre': 'loc', 'tf': 'TF'})[['TF', 'loc', 'score']]
    output_tfb.to_csv(tfb_output, sep='\t', index=False)
    
    return rna_output, peaks_output, tfb_output

def infer_grn(rna_path, peaks_path, tfb_path, grn_output):
    """Infer gene regulatory network using dictys"""
    print("Inferring gene regulatory network...")
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(grn_output), exist_ok=True)
    
    # Run dictys model to infer GRN
    os.system(f'python3 -m dictys model direct --tf_num 0 {rna_path} {peaks_path} {tfb_path} {grn_output}')
    
    return grn_output
