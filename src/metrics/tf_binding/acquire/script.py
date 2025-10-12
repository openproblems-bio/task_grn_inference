import pandas as pd
import os
import urllib.request
from tqdm import tqdm

def download_peaks_chip_atlas(celltype, genome, local_path, qval):
    def get_peak_file_name(celltype, genome, local_path, qval):
        FILE_LIST = 'fileList.tab'
        FILE_LIST_URL = f'https://chip-atlas.dbcls.jp/data/metadata/{FILE_LIST}'
        file_path = os.path.join(local_path, FILE_LIST)
        if not os.path.exists(file_path):
            os.makedirs(local_path, exist_ok=True)
            print(f'Downloading {FILE_LIST} from {FILE_LIST_URL} to {file_path}')
            urllib.request.urlretrieve(FILE_LIST_URL, file_path)
        df = pd.read_csv(file_path, sep='\t', header=None, dtype=str, comment='#', usecols=range(7))
        ct = str(celltype).strip().casefold()
        s_ct = df.iloc[:, 5].fillna('').str.strip().str.casefold()
        return df[(df[1] == genome) & (df[2] == "TFs and others") & (df[6] == qval) & (s_ct == ct)][0].iloc[0]
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
    return df

# def get_tfs_from_peaks(peaks_df, format='chip_atlas'):
#     if format == 'chip_atlas':
        
#     else:
#         raise ValueError(f"Unknown format: {format}")
    
#     return unique_tfs

def subset_peaks_by_tf(peaks_df, tf_name):
    tf_pattern = f'Name={tf_name.replace(" ", "%20")}%'
    mask = peaks_df['metadata'].str.contains(tf_pattern, na=False)
    return peaks_df[mask][['chromosome', 'start', 'end']].reset_index(drop=True)
def download_tss_file(tss_local_path):
    TSS_FILE = 'refGene.txt.gz'
    TSS_URL = f'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/{TSS_FILE}'
    TSS_PATH = os.path.join(local_path, TSS_FILE)
    if not os.path.exists(tss_local_path):
        print(f'Downloading TSS (refGene) from {TSS_URL} to {tss_local_path}')
        urllib.request.urlretrieve(TSS_URL, tss_local_path)

def genes_with_peaks_near_tss(peaks_df, tss_file, window_bp=1000):
    if peaks_df is None or len(peaks_df) == 0:
        raise ValueError("Peaks DataFrame is empty or None.")
    
    ref = pd.read_csv(
        tss_file,
        sep='\t',
        header=None,
        dtype=str,
        comment='#',
        compression='infer',
    )
    chrom = ref.iloc[:, 2].astype(str)
    strand = ref.iloc[:, 3].astype(str)
    tx_start = pd.to_numeric(ref.iloc[:, 4], errors='coerce').fillna(0).astype(int)
    tx_end = pd.to_numeric(ref.iloc[:, 5], errors='coerce').fillna(0).astype(int)
    gene = ref.iloc[:, 12].astype(str)
    tss = tx_start.where(strand == '+', tx_end)
    win_start = (tss - int(window_bp)).clip(lower=0)
    win_end = tss + int(window_bp)
    genes_df = pd.DataFrame({
        'chrom': chrom,
        'start': win_start.astype(int),
        'end': win_end.astype(int),
        'gene': gene.str.strip(),
    })
    genes_df = genes_df[(genes_df['chrom'].str.startswith('chr')) & (genes_df['gene'] != '')]
    p_chrom = peaks_df.iloc[:, 0].astype(str)
    p_start = pd.to_numeric(peaks_df.iloc[:, 1], errors='coerce').fillna(-1).astype(int)
    p_end = pd.to_numeric(peaks_df.iloc[:, 2], errors='coerce').fillna(-1).astype(int)
    peaks_simple = pd.DataFrame({'chrom': p_chrom, 'start': p_start, 'end': p_end})
    peaks_simple = peaks_simple[(peaks_simple['start'] >= 0) & (peaks_simple['end'] >= 0)]
    hits = []
    for c in genes_df['chrom'].unique():
        g = genes_df[genes_df['chrom'] == c]
        p = peaks_simple[peaks_simple['chrom'] == c]
        if p.empty or g.empty:
            continue
        for _, row in g.iterrows():
            gs, ge = row['start'], row['end']
            if ((p['end'] >= gs) & (p['start'] <= ge)).any():
                hits.append(row['gene'])
    return sorted(set(hits))

import pandas as pd
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

def _process_tf(tf, peaks_df, tss_file, window_bp):
    peaks = peaks_df[peaks_df['tf'] == tf]
    if peaks is None or len(peaks) == 0:
        return []
    genes = genes_with_peaks_near_tss(peaks, tss_file, window_bp=window_bp)
    if not genes:
        return []
    return [{'source': tf, 'target': g} for g in genes]

def build_grn(peaks_df, tss_file, genome='hg38', window_bp=1000, max_workers=None):
    tfs = peaks_df['tf'].unique()
    rows = []

    # Use ProcessPoolExecutor for parallel processing
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_process_tf, tf, peaks_df, tss_file, window_bp): tf for tf in tfs}
        for future in tqdm(as_completed(futures), total=len(futures)):
            result = future.result()
            if result:
                rows.extend(result)

    grn = pd.DataFrame(rows, columns=['source', 'target']).drop_duplicates()
    return grn
if __name__ == '__main__':
    genome = 'hg38'
    window_bp = 1000
    qval = "50"
    local_path = 'resources/datasets_raw/chipseq/'
    tss_local_path = os.path.join(local_path, 'TSS.txt.gz')
    download_tss_file(tss_local_path)
    if False:
        for cell_type in []: #'PBMC' 'K-562' 'HCT 116' 
            output_csv_path = f'resources/grn_benchmark/ground_truth/{cell_type.replace("-", "").replace(" ", "")}_chipatlas.csv'
            os.makedirs(f"{local_path}/chip_atlas/", exist_ok=True)
            print(f'Processing cell type: {cell_type}', flush=True)
            peaks_file = download_peaks_chip_atlas(cell_type, genome, f"{local_path}/chip_atlas/", qval)
            peaks_df = read_peak_file(peaks_file, source='chip_atlas')
            if True: 
                # Filter for score == 1000
                print('Filtering peaks with score == 1000', flush=True)
                print(f'Original number of peaks: {len(peaks_df)}', flush=True)
                peaks_df = peaks_df[peaks_df['score'] == '1000']
                print(f'Number of peaks after filtering: {len(peaks_df)}', flush=True)
            grn = build_grn(peaks_df, tss_local_path, genome=genome, window_bp=window_bp, max_workers=10)
            grn.to_csv(output_csv_path, index=False)
    
    if True: # remap
        if False:
            print('Processing remap2022_all_macs2_hg38_v1_0', flush=True)
            peaks_file = f"{local_path}/remap/remap2022_all_macs2_hg38_v1_0.bed.gz"
            print('Reading peaks file', flush=True)
            peaks_df = read_peak_file(peaks_file, source='remap')
            df = peaks_df.groupby('cell_type').size()
            top_cells = df.sort_values(ascending=False).head(20).index
            os.makedirs('resources/datasets_raw/chipseq/remap/', exist_ok=True)
            for cell_type in top_cells:
                subset = peaks_df[peaks_df['cell_type'] == cell_type]
                print(cell_type, subset.shape[0])
                subset.to_csv(
                    f"resources/datasets_raw/chipseq/remap/{cell_type.replace('/', '').replace('-', '').replace('.', '').replace(' ', '')}.bed",
                    sep='\t',
                    index=False,
                    header=False
                )
        # if True: 
        #     # Filter for score == 1000
        #     print('Filtering peaks with score == 1000', flush=True)
        #     print(f'Original number of peaks: {len(peaks_df)}', flush=True)
        #     peaks_df = peaks_df[peaks_df['score'] == '1000']
        #     print(f'Number of peaks after filtering: {len(peaks_df)}', flush=True)
        for cell_type in ['HEK293T', 'K-562', 'Hep-G2', 'HEK293', 'GM12878']:
            print(f'Processing cell type: {cell_type}', flush=True)
            df_sub = peaks_df[peaks_df['cell_type'] == cell_type]
            grn = build_grn(df_sub, tss_local_path, genome=genome, window_bp=window_bp)
            output_csv_path = f'resources/grn_benchmark/ground_truth/{cell_type.replace("-", "")}_remap.csv'
            print(f'Writing GRN to {output_csv_path}', flush=True)
            grn.to_csv(output_csv_path, index=False)

        