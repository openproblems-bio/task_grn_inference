"""
Build a TF→gene skeleton from TF motif data (ENCODE + JASPAR) for any dataset.

Two modes:
  - Full (ATAC available): promoter-based + peak-based TF→gene links merged
  - Motif-only (no ATAC):  promoter-based TF→gene links only (TSS ±flank_length)

Adapted from geneRNBI/src/stability_analysis/skeleton/skeleton_build.py

Usage:
    # Full mode (has ATAC):
    python build_skeleton.py --dataset op

    # Motif-only mode (no ATAC):
    python build_skeleton.py --dataset nakatake --no_atac

Outputs:
    resources/grn_benchmark/prior/skeleton_{dataset}.csv
"""

import argparse
import os
import sys
import anndata as ad
import networkx as nx
import pandas as pd
import scglue

REPO = '/home/jnourisa/projs/ongoing/task_grn_inference'
INFERENCE_DIR = f'{REPO}/resources/grn_benchmark/inference_data'
PRIOR_DIR     = f'{REPO}/resources/grn_benchmark/prior'
SUPP_DIR      = f'{REPO}/resources/supp_data'

parser = argparse.ArgumentParser()
parser.add_argument('--dataset',       type=str, required=True)
parser.add_argument('--no_atac',       action='store_true', help='Motif-only mode (no ATAC peaks)')
parser.add_argument('--flank_length',  type=int, default=1000)
parser.add_argument('--extend_range',  type=int, default=150_000)
args = parser.parse_args()

DATASET = args.dataset

par = {
    'rna':             f'{INFERENCE_DIR}/{DATASET}_rna.h5ad',
    'atac':            f'{INFERENCE_DIR}/{DATASET}_atac.h5ad',
    'annotation_file': f'{SUPP_DIR}/gencode.v47.annotation.gtf.gz',
    'motif_encode':    f'{SUPP_DIR}/databases/scglue/ENCODE-TF-ChIP-hg38.bed.gz',
    'motif_jaspar':    f'{SUPP_DIR}/databases/scglue/JASPAR2022-hg38.bed.gz',
    'temp_dir':        f'{PRIOR_DIR}/tmp/skeleton_{DATASET}',
    'out_skeleton':    f'{PRIOR_DIR}/skeleton_{DATASET}.csv',
    'flank_length':    args.flank_length,
    'extend_range':    args.extend_range,
    'no_atac':         args.no_atac,
}

os.makedirs(par['temp_dir'], exist_ok=True)


def preprocess(par):
    """Annotate RNA var with genomic coords from GTF. If ATAC, also build peak→gene guidance graph."""
    import scanpy as sc
    print('Reading RNA ...', flush=True)
    rna = ad.read_h5ad(par['rna'])

    scglue.data.get_gene_annotation(rna, gtf=par['annotation_file'], gtf_by='gene_name')
    rna = rna[:, ~rna.var['chrom'].isna()].copy()

    if 'highly_variable' not in rna.var.columns:
        sc.pp.highly_variable_genes(rna, flavor='seurat_v3', n_top_genes=min(3000, rna.n_vars))
        rna.var['highly_variable'] = rna.var['highly_variable'].fillna(False)

    str_cols = [c for c in rna.var.columns if rna.var[c].dtype == object]
    rna.var[str_cols] = rna.var[str_cols].astype(str)
    rna.write(f"{par['temp_dir']}/rna.h5ad")

    if not par['no_atac']:
        print('Reading ATAC ...', flush=True)
        atac = ad.read_h5ad(par['atac'])
        split = atac.var_names.str.split(r'[:-]')
        atac.var['chrom']      = split.map(lambda x: x[0])
        atac.var['chromStart'] = split.map(lambda x: x[1]).astype(int)
        atac.var['chromEnd']   = split.map(lambda x: x[2]).astype(int)

        guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac, extend_range=par['extend_range'])
        scglue.graph.check_graph(guidance, [rna, atac])
        atac.write(f"{par['temp_dir']}/atac.h5ad")
        nx.write_graphml(guidance, f"{par['temp_dir']}/guidance.graphml.gz")

        fwd = guidance.edge_subgraph(
            e for e, attr in dict(guidance.edges).items() if attr['type'] == 'fwd'
        ).copy()
        peak2gene = pd.DataFrame([(e[1], e[0]) for e in fwd.edges], columns=['peak', 'target'])
        peak2gene.to_csv(f"{par['temp_dir']}/peak2gene.csv", index=False)
        print(f'peak2gene: {len(peak2gene)} entries', flush=True)


def build_skeleton_promoter(par, motif_file, label):
    """TF motifs vs TSS flanks → TF→gene edges (works for all datasets)."""
    rna = ad.read_h5ad(f"{par['temp_dir']}/rna.h5ad")
    flank_bed = scglue.genomics.Bed(rna.var).strand_specific_start_site().expand(
        par['flank_length'], par['flank_length']
    )
    motif_bed = scglue.genomics.read_bed(motif_file)
    flank2tf  = scglue.genomics.window_graph(flank_bed, motif_bed, 0, right_sorted=True)
    df = pd.DataFrame([(e[1], e[0]) for e in flank2tf.edges], columns=['source', 'target'])
    out = f"{par['temp_dir']}/skeleton_promoter_{label}.csv"
    df.to_csv(out, index=False)
    print(f'Promoter skeleton ({label}): {len(df)} edges', flush=True)
    return df


def build_skeleton_peak(par, motif_file, label):
    """TF motifs vs ATAC peaks, bridged via peak→gene links."""
    atac      = ad.read_h5ad(f"{par['temp_dir']}/atac.h5ad")
    peak2gene = pd.read_csv(f"{par['temp_dir']}/peak2gene.csv")
    motif_bed = scglue.genomics.read_bed(motif_file)

    peak_bed   = scglue.genomics.Bed(atac.var)
    peak2tf    = scglue.genomics.window_graph(peak_bed, motif_bed, 0, right_sorted=True)
    df_peak2tf = pd.DataFrame([(e[1], e[0]) for e in peak2tf.edges], columns=['source', 'peak'])
    tf2gene    = peak2gene.merge(df_peak2tf, on='peak', how='inner')[['source', 'target']].drop_duplicates()
    out = f"{par['temp_dir']}/skeleton_peak_{label}.csv"
    tf2gene.to_csv(out, index=False)
    print(f'Peak skeleton ({label}): {len(tf2gene)} edges', flush=True)
    return tf2gene


if __name__ == '__main__':
    print(f'Building skeleton for {DATASET}  (no_atac={par["no_atac"]})', flush=True)

    preprocess(par)

    motif_files = {'encode': par['motif_encode'], 'jaspar': par['motif_jaspar']}

    promoter_parts = [build_skeleton_promoter(par, f, label) for label, f in motif_files.items()]
    skeleton = pd.concat(promoter_parts).drop_duplicates()

    if not par['no_atac']:
        peak_parts = [build_skeleton_peak(par, f, label) for label, f in motif_files.items()]
        skeleton_peak = pd.concat(peak_parts).drop_duplicates()
        skeleton = pd.concat([skeleton, skeleton_peak]).drop_duplicates()

    skeleton = skeleton[['source', 'target']]
    print(f'Final skeleton: {len(skeleton)} edges, '
          f'{skeleton["source"].nunique()} TFs, {skeleton["target"].nunique()} genes', flush=True)
    skeleton.to_csv(par['out_skeleton'], index=False)
    print(f'Saved to {par["out_skeleton"]}', flush=True)
