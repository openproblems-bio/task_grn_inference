import pandas as pd
import numpy as np
import anndata as ad
import os
import argparse
import warnings
warnings.filterwarnings("ignore")


parser = argparse.ArgumentParser(description="", usage="")
parser.add_argument('--rna_file', required=True)
parser.add_argument('--atac_file', required=True)
parser.add_argument('--exp_path', required=True)
parser.add_argument('--pks_path', required=True)
parser.add_argument('--frags_path', required=True)

args = vars(parser.parse_args())
rna_file = args['rna_file']
atac_file = args['atac_file']
exp_path = args['exp_path']
pks_path = args['pks_path']
frags_path = args['frags_path']

par = {
    'rna': 'resources_test/grn_benchmark/inference_data/op_rna.h5ad',
    'atac': 'resources_test/grn_benchmark/inference_data/op_atac.h5ad',
    'temp_dir': 'output/temp/'
}



# Write the RNA matrix
rna = ad.read_h5ad(rna_file)
atac = ad.read_h5ad(atac_file)
print(rna.obs.head())

def extract_exp(rna, par):
    if 'counts' not in rna.layers:
        raise ValueError('RNA data must have a "counts" layer.')
    rna_X = pd.DataFrame(np.array(rna.layers['counts'].todense()).T, columns=rna.obs.index, index=rna.var.index)
    rna_X.to_csv(par['exp_path'], sep="\t", compression="gzip")
    all_ids = rna.obs.index
    with open(par['barcodes'], "w") as f:
        for i in all_ids:
            f.write(f"{i}\n")

def extract_peaks(atac, par):
    all_atac_peak = np.unique([n.replace(':', '-') for n in atac.var.index])

    all_atac_peak = pd.DataFrame([n.split('-') for n in all_atac_peak])
    all_atac_peak.columns = ['chr', 'srt', 'end']
    all_atac_peak['srt'] = all_atac_peak['srt'].astype(int)
    all_atac_peak['end'] = all_atac_peak['end'].astype(int)
    all_atac_peak = all_atac_peak[(all_atac_peak.end - all_atac_peak.srt) >= 100]
    all_atac_peak = all_atac_peak.sort_values(by=['chr', 'srt', 'end'])
    return all_atac_peak


def process_peak(adata_atac):
    # Convert sparse matrix to COO
    coo_matrix = adata_atac.X.tocoo()
    rows = coo_matrix.row
    cols = coo_matrix.col
    counts = coo_matrix.data
    row_names = adata_atac.obs.index[rows]
    coord_names = adata_atac.var.index[cols]

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
    temp_path = frags_path.rstrip('.gz')
    df.to_csv(temp_path, sep='\t', index=False, header=False)

    # Compress TSV
    print(f'Sort and compress tsv file {frags_path}')
    os.system(f"sort -k1,1 -k2,2n {temp_path} | bgzip -c > {frags_path}")

def extract_data(par):
    par['exp_path'] = f"{par['temp_dir']}/expression.tsv.gz"
    par['barcodes'] = f"{par['temp_dir']}/barcodes_all.txt"
    par['pks_path'] = f"{par['temp_dir']}/peaks.bed"
    par['frags_path'] = f"{par['temp_dir']}/fragments.tsv.gz"
    os.makedirs(par['temp_dir'], exist_ok=True)
    
    extract_exp(rna, par)
    atac_peaks = extract_peaks(atac)
    atac_peaks.to_csv(pks_path, sep='\t', header=False, index=False)
    process_peak(atac)