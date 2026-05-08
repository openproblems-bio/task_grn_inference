"""
Build a SINGLE shared TF→gene skeleton using TF motif data (ENCODE + JASPAR).

Uses the pre-computed TSS coordinates from tss_h38.bed (protein-coding genes
only, ~19k genes) instead of a dataset-specific RNA file. This ensures the
skeleton is dataset-agnostic and covers all protein-coding genes in the human
genome.

Output: resources/grn_benchmark/prior/skeleton_motif.csv

One shared skeleton is used for ALL non-ATAC datasets (gene names and
Ensembl-ID datasets alike — it's fine if some Ensembl-ID genes are absent).

Usage:
    python skeleton_motif_only.py
    python skeleton_motif_only.py --flank_length 1000  # default
    python skeleton_motif_only.py --out_name skeleton_motif.csv
"""
import argparse, os
import anndata as ad
import pandas as pd
import numpy as np
import scglue

REPO = '/home/jnourisa/projs/ongoing/task_grn_inference'

parser = argparse.ArgumentParser()
parser.add_argument('--flank_length', type=int, default=1000)
parser.add_argument('--out_name', default='skeleton_motif.csv',
                    help='Output filename under prior/ (default: skeleton_motif.csv)')
args = parser.parse_args()

par = {
    'tss_bed':      f'{REPO}/resources/supp_data/tss_h38.bed',
    'motif_encode': f'{REPO}/resources/supp_data/databases/scglue/ENCODE-TF-ChIP-hg38.bed.gz',
    'motif_jaspar': f'{REPO}/resources/supp_data/databases/scglue/JASPAR2022-hg38.bed.gz',
    'out_skeleton': f'{REPO}/resources/grn_benchmark/prior/{args.out_name}',
    'tf_all':       f'{REPO}/resources/grn_benchmark/prior/tf_all.csv',
    'flank':        args.flank_length,
}

os.makedirs(os.path.dirname(par['out_skeleton']), exist_ok=True)

# Load known TFs for filtering sources (no header — one TF name per line)
tf_set = set(pd.read_csv(par['tf_all'], header=None, names=['TF'])['TF'].tolist())

# ── Load protein-coding TSS coordinates ─────────────────────────────────────
# tss_h38.bed has one row per TSS isoform (87k rows, 19k unique genes).
# All entries are protein_coding.
print(f"Loading TSS coordinates from {par['tss_bed']} ...", flush=True)
tss = pd.read_csv(
    par['tss_bed'], sep='\t', comment='#',
    names=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'transcript_type']
)
print(f"  {tss['name'].nunique():,} unique genes, {len(tss):,} TSS isoforms", flush=True)

# Build a Bed-compatible DataFrame: index = gene name, cols = chrom/chromStart/chromEnd/strand
# IMPORTANT: scGLUE window_graph uses the 'name' column (not the DataFrame index) as node ID.
# We must keep 'name' as both index AND column so node IDs are gene names.
gene_df = tss[['chrom', 'chromStart', 'chromEnd', 'strand', 'name']].copy()
gene_df = gene_df.set_index('name')
gene_df['name'] = gene_df.index   # keep 'name' column so window_graph uses gene names

# ── Build promoter skeleton for each motif database ─────────────────────────
parts = []
for label, motif_path in [('encode', par['motif_encode']), ('jaspar', par['motif_jaspar'])]:
    print(f"Building promoter skeleton ({label}) ...", flush=True)
    motif_bed = scglue.genomics.read_bed(motif_path)

    # strand_specific_start_site collapses each gene to its TSS (respects strand),
    # then expand by flank_length on both sides
    flank_bed = scglue.genomics.Bed(gene_df).strand_specific_start_site().expand(par['flank'], par['flank'])

    flank2tf = scglue.genomics.window_graph(flank_bed, motif_bed, 0, right_sorted=True)

    df = pd.DataFrame(
        [(e[1], e[0]) for e in flank2tf.edges],
        columns=['source', 'target']
    )
    # Filter sources to known TFs only
    df = df[df['source'].isin(tf_set)]
    print(f"  {label}: {len(df):,} TF→gene edges (before dedup)", flush=True)
    parts.append(df)

skeleton = pd.concat(parts).drop_duplicates()
print(
    f"\nFinal skeleton: {len(skeleton):,} edges, "
    f"{skeleton['source'].nunique():,} TFs, {skeleton['target'].nunique():,} genes",
    flush=True
)
skeleton.to_csv(par['out_skeleton'], index=False)
print(f"Saved: {par['out_skeleton']}", flush=True)
