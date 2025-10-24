import pandas as pd
import os
import numpy as np
import urllib.request
from tqdm import tqdm

# panda prints all the columns
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 1000)

def download_peaks_chip_atlas(celltype, genome, local_path, qval):
    """Download ChIP-Atlas peaks for a specific cell type and genome."""
    def get_peak_file_name(celltype, genome, local_path, qval):
        FILE_LIST = 'fileList.tab'
        FILE_LIST_URL = f'https://chip-atlas.dbcls.jp/data/metadata/{FILE_LIST}'
        file_path = os.path.join(local_path, FILE_LIST)
        if not os.path.exists(file_path):
            os.makedirs(local_path, exist_ok=True)
            print(f'Downloading {FILE_LIST} from {FILE_LIST_URL} to {file_path}')
            urllib.request.urlretrieve(FILE_LIST_URL, file_path)
        df = pd.read_csv(file_path, sep='\t', header=None, dtype=str, comment='#', usecols=range(7))
        assert df.shape[0] > 0, f'No entries found in {file_path}'
        ct = str(celltype).strip().casefold()
        s_ct = df.iloc[:, 5].fillna('').str.strip().str.casefold()
        
        df_sub = df[(df[1] == genome) & (df[2] == "TFs and others") & (df[6] == qval) & (s_ct == ct)][0]
        
        if df_sub.shape[0] == 0:
            raise ValueError(f'No peaks found for cell type {celltype}, {genome}, {qval} in {file_path}')
        return df_sub.iloc[0]
    
    filename = get_peak_file_name(celltype, genome, local_path, qval)
    peaks_path = os.path.join(local_path, filename + '.bed')
    if not os.path.exists(peaks_path):
        os.makedirs(local_path, exist_ok=True)
        PEAKS_URL = f'https://chip-atlas.dbcls.jp/data/{genome}/assembled/'
        peaks_url = PEAKS_URL + filename + '.bed'
        print(f'Downloading peaks from {peaks_url} to {peaks_path}')
        urllib.request.urlretrieve(peaks_url, peaks_path)
    return peaks_path

def read_peak_file(peaks_path, source):
    """Read and process peak files from different sources."""
    df = pd.read_csv(
        peaks_path,
        sep='\t',
        header=None,
        dtype=str,
        comment='#',
        names=['chromosome', 'start', 'end', 'metadata', 'score', 'strand', 
               'thick_start', 'thick_end', 'rgb']
    )
    if source == 'chip_atlas':
        tf_names = df['metadata'].str.extract(r'Name=([^%]+)%')[0]
        tfs = (tf_names
                 .dropna()
                 .str.strip()
                 .str.replace('%20', ' '))
        df['tf'] = tfs
    elif source == 'remap':
        df[['experiment_id', 'tf', 'cell_type']] = df['metadata'].str.split('.', n=2, expand=True)
    else:
        raise ValueError(f"Unknown source: {source}")
    df['score'] = pd.to_numeric(df['score'], errors='coerce').fillna(0).astype(float)
    return df

def subset_peaks_by_tf(peaks_df, tf_name):
    """Subset peaks DataFrame by transcription factor name."""
    tf_pattern = f'Name={tf_name.replace(" ", "%20")}%'
    mask = peaks_df['metadata'].str.contains(tf_pattern, na=False)
    return peaks_df[mask][['chromosome', 'start', 'end']].reset_index(drop=True)

def download_tss_file(tss_local_path):
    """Download TSS (Transcription Start Sites) file from UCSC."""
    TSS_FILE = 'refGene.txt.gz'
    TSS_URL = f'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/{TSS_FILE}'

    if not os.path.exists(tss_local_path):
        print(f'Downloading TSS (refGene) from {TSS_URL} to {tss_local_path}')
        os.makedirs(os.path.dirname(tss_local_path), exist_ok=True)
        urllib.request.urlretrieve(TSS_URL, tss_local_path)

def download_and_prepare_chip_atlas_data(cell_types, genome='hg38', qval="50", local_path='resources/datasets_raw/chipseq/'):
    """Download and prepare ChIP-Atlas data for multiple cell types."""
    print('Processing chip-atlas', flush=True)
    prepared_data = {}
    
    for cell_type in cell_types:
        cell_type_c = cell_type.replace('/', '').replace('-', '').replace('.', '').replace(' ', '').replace('_', '').upper()
        if cell_type_c == 'HEK293TREx':
            cell_type_c = 'HEK293T'
        
        os.makedirs(f"{local_path}/chip_atlas/", exist_ok=True)
        print(f'Processing cell type: {cell_type}', flush=True)
        
        peaks_file = download_peaks_chip_atlas(cell_type, genome, f"{local_path}/chip_atlas/", qval)
        peaks_df = read_peak_file(peaks_file, source='chip_atlas')
        
        # Filter for score > 50
        print('Filtering peaks with score > 50', flush=True)
        print(f'Original number of peaks: {len(peaks_df)}', flush=True)
        peaks_df = peaks_df[peaks_df['score'] > 50]
        print(f'Number of peaks after filtering: {len(peaks_df)}', flush=True)
        
        prepared_data[cell_type_c] = peaks_df
        
    return prepared_data

def download_and_prepare_remap_data(cell_types, local_path='resources/datasets_raw/chipseq/'):
    """Download and prepare ReMap data for multiple cell types."""
    print('Processing remap2', flush=True)
    prepared_data = {}
    
    # First, process the main remap file to extract cell-type specific data
    print('Processing remap2022_all_macs2_hg38_v1_0', flush=True)
    peaks_file = f"{local_path}/remap/remap2022_all_macs2_hg38_v1_0.bed.gz"
    print('Reading peaks file', flush=True)
    peaks_df = read_peak_file(peaks_file, source='remap')
    df_grouped = peaks_df.groupby('cell_type').size()
    top_cells = df_grouped.sort_values(ascending=False).head(20).index
    
    os.makedirs('resources/datasets_raw/chipseq/remap/', exist_ok=True)
    for cell_type in top_cells:
        subset = peaks_df[peaks_df['cell_type'] == cell_type]
        print(cell_type, subset.shape[0])
        
        subset.to_csv(
            f"resources/datasets_raw/chipseq/remap/{cell_type}.bed",
            sep='\t',
            index=False,
            header=False
        )
    
    # Now process the specific cell types
    for cell_type in cell_types:
        if cell_type not in top_cells:
            print(f'Warning: {cell_type} not found in top cell types')
            continue
            
        print(f'Processing cell type: {cell_type}', flush=True)
        peaks_file = f"resources/datasets_raw/chipseq/remap/{cell_type}.bed"
        df = pd.read_csv(
            peaks_file,
            sep='\t',
            header=None,
            dtype=str,
            comment='#',
            names=['chromosome', 'start', 'end', 'metadata', 'score', 'strand', 
                'thick_start', 'thick_end', 'rgb', '_', 'tf', 'cell_type']
        )
        df['score'] = pd.to_numeric(df['score'], errors='coerce')

        print('Shape before filtering:', df.shape)
        print("Before filtering, ", df.shape[0], " peaks")
        df = df[df['score'] > 0]
        print("After filtering, ", df.shape[0], " peaks")
        print('Shape after filtering:', df.shape)
        
        prepared_data[cell_type] = df
        
    return prepared_data

def main():
    """Main function to download and prepare data."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Download and prepare ChIP-seq data')
    parser.add_argument('--genome', default='hg38', help='Genome version (default: hg38)')
    parser.add_argument('--qval', default='50', help='Q-value threshold (default: 50)')
    parser.add_argument('--local-path', default='resources/datasets_raw/chipseq/', 
                       help='Local path for data storage')
    parser.add_argument('--chip-atlas-cells', nargs='+', default=['k-562'], 
                       help='ChIP-Atlas cell types to process')
    parser.add_argument('--remap-cells', nargs='+', default=['K562'], 
                       help='ReMap cell types to process')
    parser.add_argument('--skip-chip-atlas', action='store_true', 
                       help='Skip ChIP-Atlas data processing')
    parser.add_argument('--skip-remap', action='store_true', 
                       help='Skip ReMap data processing')
    
    args = parser.parse_args()
    
    tss_local_path = os.path.join(args.local_path, 'TSS.txt.gz')
    
    print("=== Data Download and Preparation Pipeline ===")
    
    # Download TSS file
    print("\n1. Downloading TSS file...")
    download_tss_file(tss_local_path)
    
    chip_atlas_data = {}
    remap_data = {}
    
    # Download ChIP-Atlas data
    if not args.skip_chip_atlas and args.chip_atlas_cells:
        print("\n2. Processing ChIP-Atlas data...")
        chip_atlas_data = download_and_prepare_chip_atlas_data(
            args.chip_atlas_cells, args.genome, args.qval, args.local_path
        )
    
    # Download ReMap data
    if not args.skip_remap and args.remap_cells:
        print("\n3. Processing ReMap data...")
        remap_data = download_and_prepare_remap_data(args.remap_cells, args.local_path)
    
    print("\n=== Data download and preparation completed! ===")
    if chip_atlas_data:
        print(f"ChIP-Atlas data prepared for: {list(chip_atlas_data.keys())}")
    if remap_data:
        print(f"ReMap data prepared for: {list(remap_data.keys())}")
    
    return chip_atlas_data, remap_data

if __name__ == '__main__':
    main()