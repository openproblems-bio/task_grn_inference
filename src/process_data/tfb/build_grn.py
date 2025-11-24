import pandas as pd
import os
import numpy as np
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

# panda prints all the columns
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 1000)

def genes_with_peaks_near_tss(peaks_df, tss_file, window_up=1000, window_down=100):
    """Find genes that have peaks near their transcription start sites."""
    if peaks_df is None or len(peaks_df) == 0:
        raise ValueError("Peaks DataFrame is empty or None.")
    
    # Check if tss_file is already a DataFrame or a file path
    if isinstance(tss_file, str):
        # Check if it's a BED format file (has header starting with #)
        with open(tss_file, 'r') as f:
            first_line = f.readline().strip()
        
        if first_line.startswith('#'):
            # BED format - manually skip comment lines because pandas comment handling is buggy
            lines = []
            with open(tss_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        lines.append(line.strip())
            
            # Parse lines into dataframe
            data = []
            for line in lines:
                data.append(line.split('\t'))
            
            # Create DataFrame with proper column names
            ref = pd.DataFrame(data, columns=['Chromosome', 'Start', 'End', 'Gene', 'Score', 'Strand', 'Transcript_type'])
            
            chrom = ref['Chromosome'].astype(str)
            strand = ref['Strand'].astype(str) 
            tx_start = pd.to_numeric(ref['Start'], errors='coerce').fillna(0).astype(int)
            tx_end = pd.to_numeric(ref['End'], errors='coerce').fillna(0).astype(int)
            gene = ref['Gene'].astype(str)
        else:
            # Original format (refGene)
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
    else:
        # DataFrame passed directly
        ref = tss_file
        chrom = ref.iloc[:, 2].astype(str)
        strand = ref.iloc[:, 3].astype(str)
        tx_start = pd.to_numeric(ref.iloc[:, 4], errors='coerce').fillna(0).astype(int)
        tx_end = pd.to_numeric(ref.iloc[:, 5], errors='coerce').fillna(0).astype(int)
        gene = ref.iloc[:, 12].astype(str)
    tss = tx_start.where(strand == '+', tx_end)
    win_start = np.where(strand == '+', tss - window_up, tss - window_down)
    win_end = np.where(strand == '+', tss + window_down, tss + window_up)
    genes_df = pd.DataFrame({
        'chrom': chrom,
        'start': np.clip(win_start, 0, None).astype(int), 
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

def _process_tf(tf, peaks_df, tss_file):
    """Process a single transcription factor to find target genes."""
    peaks = peaks_df[peaks_df['tf'] == tf]
    if peaks is None or len(peaks) == 0:
        return []
    genes = genes_with_peaks_near_tss(peaks, tss_file)
    if not genes:
        return []
    return [{'source': tf, 'target': g} for g in genes]

def build_grn(peaks_df, tss_file, genome='hg38', num_workers=None):
    """Build Gene Regulatory Network from peaks and TSS data."""
    tfs = peaks_df['tf'].unique()
    rows = []

    # Use ProcessPoolExecutor for parallel processing
    with ProcessPoolExecutor(num_workers=num_workers) as executor:
        futures = {executor.submit(_process_tf, tf, peaks_df, tss_file): tf for tf in tfs}
        for future in tqdm(as_completed(futures), total=len(futures)):
            result = future.result()
            if result:
                rows.extend(result)

    grn = pd.DataFrame(rows, columns=['source', 'target']).drop_duplicates()
    return grn


def read_tfb_file(tfb_file, source_type):
    """Read TF binding files and convert to the format expected by build_grn functions."""
    # Read file without setting column names first
    df = pd.read_csv(tfb_file, sep='\t', header=None, comment='#')
    
    print(f"File has {len(df.columns)} columns")
    print(f"First row: {df.iloc[0].tolist()}")
    
    # Assign column names based on actual number of columns and source type
    if source_type == 'unibind':
        # Flexible column assignment for unibind
        if len(df.columns) >= 4:
            df = df.rename(columns={0: 'chromosome', 1: 'start', 2: 'end', 3: 'tf'})
            # Additional columns are optional (cell_type, etc.)
        else:
            raise ValueError(f"Expected at least 4 columns for unibind format, got {len(df.columns)}")
    elif source_type in ['chipatlas', 'remap2022']:
        # Flexible column assignment for chipatlas and remap2022
        if len(df.columns) >= 4:
            df = df.rename(columns={0: 'chromosome', 1: 'start', 2: 'end', 3: 'tf'})
            # Additional columns are optional
        else:
            raise ValueError(f"Expected at least 4 columns for {source_type} format, got {len(df.columns)}")
    else:
        raise ValueError(f"Unknown source type: {source_type}")
    
    # Convert coordinates to numeric
    df['start'] = pd.to_numeric(df['start'], errors='coerce')
    df['end'] = pd.to_numeric(df['end'], errors='coerce')
    
    # Add required columns for compatibility with existing functions
    df['score'] = 1.0  # Default score
    df['strand'] = '.'
    df['thick_start'] = df['start']
    df['thick_end'] = df['end'] 
    df['rgb'] = '0,0,0'
    df['metadata'] = df['tf']  # Use TF name as metadata
    
    return df

def read_tss_file_bed(tss_file):
    """Read TSS file in BED format and convert to format expected by genes_with_peaks_near_tss."""
    df = pd.read_csv(tss_file, sep='\t', comment='#', header=0)
    
    # Convert to the format expected by the existing function
    # The function expects columns at specific indices
    result_df = pd.DataFrame()
    result_df[0] = ''  # bin (not used)
    result_df[1] = ''  # name (not used) 
    result_df[2] = df.iloc[:, 0]  # chrom
    result_df[3] = df.iloc[:, 5]  # strand
    result_df[4] = df.iloc[:, 1]  # txStart
    result_df[5] = df.iloc[:, 2]  # txEnd
    for i in range(6, 12):
        result_df[i] = ''  # other columns not used
    result_df[12] = df.iloc[:, 3]  # gene name
    
    return result_df

def print_grn_summary(grn_results, data_source):
    """Print summary statistics for GRN results."""
    print(f"\n=== {data_source} GRN Summary ===")
    for cell_type, result in grn_results.items():
        print(f"{cell_type}:")
        print(f"  - Edges: {result['num_edges']}")
        print(f"  - TFs: {result['num_tfs']}")
        print(f"  - Target genes: {result['num_targets']}")
        print(f"  - Output: {result['output_path']}")

def main():
    """Main function to build GRNs from existing prepared data files."""
    import argparse
    import pickle
    
    parser = argparse.ArgumentParser(description='Build Gene Regulatory Networks from prepared data')
    parser.add_argument('--tss-file', required=True, help='Path to TSS file')
    parser.add_argument('--chip-atlas-data', help='Path to pickled ChIP-Atlas prepared data')
    parser.add_argument('--remap-data', help='Path to pickled ReMap prepared data')
    parser.add_argument('--peaks-file', help='Path to individual peaks file')
    parser.add_argument('--source', choices=['chip_atlas', 'remap'], 
                       help='Source type for individual peaks file')
    parser.add_argument('--output-dir', default='resources/grn_benchmark/ground_truth/', 
                       help='Output directory for GRN files')
    parser.add_argument('--max-workers', type=int, default=10, 
                       help='Maximum number of parallel workers')
    parser.add_argument('--cell-type-name', help='Cell type name for individual file processing')
    
    args = parser.parse_args()
    
    print("=== GRN Building Pipeline ===")
    
    if not os.path.exists(args.tss_file):
        print(f"Error: TSS file not found: {args.tss_file}")
        print("Please run download_data.py first or provide a valid TSS file path.")
        return
    
    results = {}
    
    # Process ChIP-Atlas data
    if args.chip_atlas_data:
        print(f"\n1. Loading ChIP-Atlas data from {args.chip_atlas_data}...")
        with open(args.chip_atlas_data, 'rb') as f:
            chip_atlas_data = pickle.load(f)
        
        print("2. Building GRNs from ChIP-Atlas data...")
        chip_atlas_results = build_grn_from_chip_atlas_data(
            chip_atlas_data, args.tss_file, args.output_dir, args.num_workers
        )
        print_grn_summary(chip_atlas_results, "ChIP-Atlas")
        results.update(chip_atlas_results)
    
    # Process ReMap data
    if args.remap_data:
        print(f"\n3. Loading ReMap data from {args.remap_data}...")
        with open(args.remap_data, 'rb') as f:
            remap_data = pickle.load(f)
        
        print("4. Building GRNs from ReMap data...")
        remap_results = build_grn_from_remap_data(
            remap_data, args.tss_file, args.output_dir, args.num_workers
        )
        print_grn_summary(remap_results, "ReMap")
        results.update(remap_results)
    
    # Process individual peaks file
    if args.peaks_file and args.source and args.cell_type_name:
        print(f"\n5. Processing individual peaks file: {args.peaks_file}...")
        from download_data import read_peak_file
        
        peaks_df = read_peak_file(args.peaks_file, args.source)
        prepared_data = {args.cell_type_name: peaks_df}
        
        if args.source == 'chip_atlas':
            individual_results = build_grn_from_chip_atlas_data(
                prepared_data, args.tss_file, args.output_dir, args.num_workers
            )
        else:  # remap
            individual_results = build_grn_from_remap_data(
                prepared_data, args.tss_file, args.output_dir, args.num_workers
            )
        
        print_grn_summary(individual_results, f"Individual {args.source}")
        results.update(individual_results)
    
    if not results:
        print("\nNo data processed. Please provide:")
        print("  - ChIP-Atlas data file (--chip-atlas-data)")
        print("  - ReMap data file (--remap-data)")
        print("  - Individual peaks file (--peaks-file + --source + --cell-type-name)")
        print("\nExample usage:")
        print("  python build_grn.py --tss-file data/TSS.txt.gz --peaks-file peaks.bed --source chip_atlas --cell-type-name K562")
        return
    
    print(f"\n=== GRN building completed! Built {len(results)} GRNs ===")

if __name__ == '__main__':
    main()