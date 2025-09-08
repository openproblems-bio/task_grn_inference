import pandas as pd
import numpy as np

import os
import argparse
import warnings
warnings.filterwarnings("ignore")

def define_vars(par):
    par['exp_path'] = f"{par['temp_dir']}/expression.tsv.gz"
    par['barcodes'] = f"{par['temp_dir']}/barcodes_all.txt"
    par['pks_path'] = f"{par['temp_dir']}/peaks.bed"
    par['frags_path'] = f"{par['temp_dir']}/fragments.tsv.gz"

    par['bam_name'] = f"{par['temp_dir']}/reads_all.bam"
    par['bai_name'] = f"{par['temp_dir']}/reads_all.bai"
    os.makedirs(par['temp_dir'], exist_ok=True)

def extract_exp(par):
    print('Extracting expression matrix and barcodes', flush=True)
    import anndata as ad
    rna = ad.read_h5ad(par['rna'])
    if 'counts' not in rna.layers:
        X = rna.X
    else:
        X = rna.layers['counts']
    if not isinstance(X, np.ndarray):
        X = X.todense()
    rna_X = pd.DataFrame(X.T, columns=rna.obs.index, index=rna.var.index)
    rna_X.to_csv(par['exp_path'], sep="\t", compression="gzip")
    all_ids = rna.obs.index
    with open(par['barcodes'], "w") as f:
        for i in all_ids:
            f.write(f"{i}\n")

def extract_atac(par):
    print('Extracting peaks and fragments', flush=True)
    import anndata as ad
    atac = ad.read_h5ad(par['atac'])
    all_atac_peak = np.unique([n.replace(':', '-') for n in atac.var.index])
    all_atac_peak = pd.DataFrame([n.split('-') for n in all_atac_peak])
    all_atac_peak.columns = ['chr', 'srt', 'end']
    all_atac_peak['srt'] = all_atac_peak['srt'].astype(int)
    all_atac_peak['end'] = all_atac_peak['end'].astype(int)
    all_atac_peak = all_atac_peak[(all_atac_peak.end - all_atac_peak.srt) >= 100]
    all_atac_peak = all_atac_peak.sort_values(by=['chr', 'srt', 'end'])
    all_atac_peak.to_csv(par['pks_path'], sep='\t', header=False, index=False)


    # get fragments
    coo_matrix = atac.X.tocoo()
    rows = coo_matrix.row
    cols = coo_matrix.col
    counts = coo_matrix.data
    row_names = atac.obs.index[rows]
    coord_names = atac.var.index[cols]

    # Build DataFrame
    d = {
        'chromosome': [s.split(':')[0] for s in coord_names],
        'start': [int(s.split(':')[1].split('-')[0]) for s in coord_names],
        'end': [int(s.split(':')[1].split('-')[1]) for s in coord_names],
        # 'barcode': [barcode.replace('-', '') for barcode in row_names],
        'barcode': [barcode for barcode in row_names],
        'count': np.squeeze(np.asarray(counts.astype(np.uint16)))
    }
    df = pd.DataFrame(d, copy=False)

    # Write to TSV without header/index
    frags_path = par['frags_path']
    temp_path = frags_path.rstrip('.gz')
    df.to_csv(temp_path, sep='\t', index=False, header=False)

    # Compress TSV
    print(f'Sort and compress tsv file {frags_path}')
    os.system(f"sort -k1,1 -k2,2n {temp_path} | bgzip -c > {frags_path}")

def create_bam(par):
    print('Creating BAM file from fragments', flush=True)
    import subprocess
    print("Creating BAM")
    cmd = f"python src/methods/dictys/frag_to_bam.py --fnames {par['frags_path']} --barcodes {par['barcodes']}"
    view_sort_cmd = f"samtools view -b | samtools sort -o {par['bam_name']}"
    full_cmd = f"{cmd} | {view_sort_cmd}"
    subprocess.run(full_cmd, shell=True, check=True)
    subprocess.run(["samtools", "index", par['bam_name'], par['bai_name']], check=True)
    
def format_inputs(par):
    # extract_exp(par)
    # extract_atac(par)
    create_bam(par)

def export_net():
    d0=network.from_file('output/static.h5') # check how to control this
    d0.export('net',sparsities=[None])  # no way to control the name it's exporting
    net = pd.read_csv('output/dictys_test/net/Full/Subset1.tsv.gz', sep='\t') # this is customized

    prediction = net.reset_index().melt(id_vars=net.index.name or 'index', 
                                 var_name='target', 
                                 value_name='weight')
    prediction.rename(columns={net.index.name or 'index': 'source'}, inplace=True)

    prediction = process_links(prediction)

    # put this into grn format