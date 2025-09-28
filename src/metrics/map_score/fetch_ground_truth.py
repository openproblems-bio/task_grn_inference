import pandas as pd
import os
import urllib.request

TAB_PATH = './data/chip_atlas/experimentList.tab'
TAB_URL = 'https://chip-atlas.dbcls.jp/data/metadata/experimentList.tab'
PEAKS_PATH = './data/chip_atlas/allPeaks_light.hg38.05.bed.gz'
PEAKS_URL = 'https://chip-atlas.dbcls.jp/data/hg38/allPeaks_light/allPeaks_light.hg38.05.bed.gz'
TSS_PATH = './data/chip_atlas/refGene.hg38.txt.gz'
TSS_URL = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz'
GENOME = 'hg38'

def _load_experiment_tab_df():
    if not os.path.exists(TAB_PATH):
        os.makedirs(os.path.dirname(TAB_PATH), exist_ok=True)
        print(f'Downloading experimentList.tab from {TAB_URL} to {TAB_PATH}')
        urllib.request.urlretrieve(TAB_URL, TAB_PATH)
    df = pd.read_csv(TAB_PATH, sep='\t', header=None, dtype=str, comment='#', usecols=range(9))
    return df[df[2] == "TFs and others"]


def _mask_celltype_and_genome(df, cell_type):
    ct = str(cell_type).strip().casefold()
    genome = str(GENOME).strip().casefold()
    s_ct = df.iloc[:, 5].fillna('').str.strip().str.casefold()
    s_gen = df.iloc[:, 1].fillna('').str.strip().str.casefold()
    return (s_ct == ct) & (s_gen == genome)


def find_experiment_ids_by_cell_and_tf(cell_type, tf_name):
    df = _load_experiment_tab_df()
    tf = str(tf_name).strip().casefold()
    s_tf = df.iloc[:, 3].fillna('').str.strip().str.casefold()
    mask = _mask_celltype_and_genome(df, cell_type)
    return (
        df.loc[mask & (s_tf == tf), 0]
        .dropna()
        .astype(str)
        .str.strip()
        .drop_duplicates()
        .tolist()
    )

def subset_peaks_by_experiment_ids(experiment_ids):
    if not os.path.exists(PEAKS_PATH):
        os.makedirs(os.path.dirname(PEAKS_PATH), exist_ok=True)
        print(f'Downloading peaks from {PEAKS_URL} to {PEAKS_PATH}')
        urllib.request.urlretrieve(PEAKS_URL, PEAKS_PATH)
    if not experiment_ids:
        return pd.DataFrame(columns=[0, 1, 2, 3, 4])
    ids = {str(x).strip() for x in experiment_ids}
    df = pd.read_csv(
        PEAKS_PATH,
        sep='\t',
        header=None,
        dtype=str,
        comment='#',
        compression='infer',
    )
    srx = df.iloc[:, 3].fillna('').str.strip()
    return df.loc[srx.isin(ids)].reset_index(drop=True)


def genes_with_peaks_near_tss(peaks_df, window_bp=1000):
    if peaks_df is None or len(peaks_df) == 0:
        return []
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


def list_tfs_for_cell_type(cell_type):
    df = _load_experiment_tab_df()
    mask = _mask_celltype_and_genome(df, cell_type)
    return (
        df.loc[mask, 3]
        .dropna()
        .astype(str)
        .str.strip()
        .drop_duplicates()
        .sort_values()
        .tolist()
    )


def build_celltype_grn(cell_type, window_bp=1000):
    safe_ct = str(cell_type).strip().replace(' ', '_')
    output_csv_path = f'./data/chip_atlas/grn_{safe_ct}_{GENOME}.csv'
    
    if os.path.exists(output_csv_path):
        return output_csv_path
        
    tfs = list_tfs_for_cell_type(cell_type)
    rows = []
    for tf in tfs:
        exp_ids = find_experiment_ids_by_cell_and_tf(cell_type, tf)
        if not exp_ids:
            continue
        peaks = subset_peaks_by_experiment_ids(exp_ids)
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

