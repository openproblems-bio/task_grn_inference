"""
Build a TF→gene skeleton using only TF motif data (ENCODE + JASPAR).
No ATAC data needed — uses promoter/TSS flanking regions only.

Output: resources/grn_benchmark/prior/skeleton_{dataset}.csv

Usage:
    python skeleton_motif_only.py --dataset nakatake
    python skeleton_motif_only.py --dataset xaira_HEK293T  # Ensembl IDs → gtf_by=gene_id
"""
import argparse, os, sys
import anndata as ad
import pandas as pd
import scglue

REPO = '/home/jnourisa/projs/ongoing/task_grn_inference'

DATASET_GTF_BY = {
    'xaira_HEK293T': 'gene_id',
    'xaira_HCT116':  'gene_id',
}

parser = argparse.ArgumentParser()
parser.add_argument('--dataset', required=True)
parser.add_argument('--flank_length', type=int, default=1000)
args = parser.parse_args()

DATASET = args.dataset
gtf_by  = DATASET_GTF_BY.get(DATASET, 'gene_name')

par = {
    'rna':          f'{REPO}/resources/grn_benchmark/inference_data/{DATASET}_rna.h5ad',
    'annotation':   f'{REPO}/resources/supp_data/gencode.v47.annotation.gtf.gz',
    'motif_encode': f'{REPO}/resources/supp_data/databases/scglue/ENCODE-TF-ChIP-hg38.bed.gz',
    'motif_jaspar': f'{REPO}/resources/supp_data/databases/scglue/JASPAR2022-hg38.bed.gz',
    'out_skeleton': f'{REPO}/resources/grn_benchmark/prior/skeleton_{DATASET}.csv',
    'tf_all':       f'{REPO}/resources/grn_benchmark/prior/tf_all.csv',
    'flank':        args.flank_length,
}

os.makedirs(os.path.dirname(par['out_skeleton']), exist_ok=True)

print(f"Loading RNA: {par['rna']}", flush=True)
rna = ad.read_h5ad(par['rna'])
print(f"  Shape: {rna.shape}", flush=True)

print(f"Annotating genes from GTF (by={gtf_by}) ...", flush=True)
scglue.data.get_gene_annotation(rna, gtf=par['annotation'], gtf_by=gtf_by)
rna = rna[:, ~rna.var['chrom'].isna()].copy()
print(f"  Genes with coords: {rna.n_vars}", flush=True)

# Load known TFs for filtering
tf_all = pd.read_csv(par['tf_all'])['TF'].tolist()
tf_set = set(tf_all)

parts = []
for label, motif_path in [('encode', par['motif_encode']), ('jaspar', par['motif_jaspar'])]:
    print(f"Building promoter skeleton ({label}) ...", flush=True)
    motif_bed = scglue.genomics.read_bed(motif_path)

    flank_bed = (scglue.genomics.Bed(rna.var)
                 .strand_specific_start_site()
                 .expand(par['flank'], par['flank']))

    flank2tf = scglue.genomics.window_graph(flank_bed, motif_bed, 0, right_sorted=True)

    df = pd.DataFrame(
        [(e[1], e[0]) for e in flank2tf.edges if e[1] in tf_set],
        columns=['source', 'target']
    )
    print(f"  {label}: {len(df)} TF→gene edges", flush=True)
    parts.append(df)

skeleton = pd.concat(parts).drop_duplicates()
print(f"Final skeleton: {len(skeleton)} edges, "
      f"{skeleton['source'].nunique()} TFs, {skeleton['target'].nunique()} genes", flush=True)
skeleton.to_csv(par['out_skeleton'], index=False)
print(f"Saved: {par['out_skeleton']}", flush=True)
