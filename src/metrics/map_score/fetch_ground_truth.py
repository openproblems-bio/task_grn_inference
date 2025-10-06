import pandas as pd
import os
import urllib.request

GENOME = 'hg38'
LOCAL_PATH = './data/chip_atlas/'
FILE_LIST = 'fileList.tab'
FILE_LIST_URL = f'https://chip-atlas.dbcls.jp/data/metadata/{FILE_LIST}'
PEAKS_URL = f'https://chip-atlas.dbcls.jp/data/{GENOME}/assembled/'
TSS_FILE = 'refGene.txt.gz'
TSS_URL = f'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/{TSS_FILE}'

def get_peak_file_name(celltype, qval="50"):
    file_path = os.path.join(LOCAL_PATH, FILE_LIST)
    if not os.path.exists(file_path):
        os.makedirs(LOCAL_PATH, exist_ok=True)
        print(f'Downloading {FILE_LIST} from {FILE_LIST_URL} to {file_path}')
        urllib.request.urlretrieve(FILE_LIST_URL, file_path)
    df = pd.read_csv(file_path, sep='\t', header=None, dtype=str, comment='#', usecols=range(7))
    ct = str(celltype).strip().casefold()
    s_ct = df.iloc[:, 5].fillna('').str.strip().str.casefold()
    return df[(df[1] == GENOME) & (df[2] == "TFs and others") & (df[6] == qval) & (s_ct == ct)][0].iloc[0]

def download_peaks_for_celltype(celltype):
    filename = get_peak_file_name(celltype)
    peaks_path = os.path.join(LOCAL_PATH, filename + '.bed')
    if not os.path.exists(peaks_path):
        os.makedirs(LOCAL_PATH, exist_ok=True)
        peaks_url = PEAKS_URL + filename + '.bed'
        print(f'Downloading peaks from {peaks_url} to {peaks_path}')
        urllib.request.urlretrieve(peaks_url, peaks_path)
    df = pd.read_csv(
        peaks_path,
        sep='\t',
        header=None,
        dtype=str,
        comment='#',
        names=['chromosome', 'start', 'end', 'metadata', 'score', 'strand', 
               'thick_start', 'thick_end', 'rgb']
    )
    return df

def get_tfs_from_peaks(peaks_df):
    tf_names = peaks_df['metadata'].str.extract(r'Name=([^%]+)%')[0]
    unique_tfs = (tf_names
                 .dropna()
                 .str.strip()
                 .str.replace('%20', ' ')
                 .drop_duplicates()
                 .sort_values()
                 .tolist())
    return unique_tfs

def subset_peaks_by_tf(peaks_df, tf_name):
    tf_pattern = f'Name={tf_name.replace(" ", "%20")}%'
    mask = peaks_df['metadata'].str.contains(tf_pattern, na=False)
    return peaks_df[mask][['chromosome', 'start', 'end']].reset_index(drop=True)

def genes_with_peaks_near_tss(peaks_df, window_bp=1000):
    if peaks_df is None or len(peaks_df) == 0:
        return []
    TSS_PATH = os.path.join(LOCAL_PATH, TSS_FILE)
    if not os.path.exists(TSS_PATH):
        os.makedirs(os.path.dirname(TSS_PATH), exist_ok=True)
        print(f'Downloading TSS (refGene) from {TSS_URL} to {TSS_PATH}')
        urllib.request.urlretrieve(TSS_URL, TSS_PATH)
    ref = pd.read_csv(
        TSS_PATH,
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

def build_celltype_grn(cell_type, window_bp=1000):
    safe_ct = str(cell_type).strip().replace(' ', '_')
    output_csv_path = LOCAL_PATH + f'grn_{safe_ct}_{GENOME}.csv'
    if os.path.exists(output_csv_path):
        return output_csv_path
    peaks_df = download_peaks_for_celltype(cell_type)
    tfs = get_tfs_from_peaks(peaks_df)
    rows = []
    for tf in tfs:
        peaks = subset_peaks_by_tf(peaks_df, tf)
        if peaks is None or len(peaks) == 0:
            continue
        genes = genes_with_peaks_near_tss(peaks, window_bp=window_bp)
        if not genes:
            continue
        rows.extend({'source': tf, 'target': g} for g in genes)
    grn = pd.DataFrame(rows, columns=['source', 'target']).drop_duplicates()
    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
    grn.to_csv(output_csv_path, index=False)
    return output_csv_path
