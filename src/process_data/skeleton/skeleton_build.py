"""
Build a TF→gene skeleton from TF motif data (ENCODE + JASPAR) using scGLUE.
Requires {dataset}_rna.h5ad and {dataset}_atac.h5ad in $INFERENCE_DIR.
Output: $RESULTS_DIR/experiment/skeleton/{dataset}/skeleton.csv

Two complementary strategies are merged:
  1. Promoter-based: TF motifs vs TSS flanks (±1 kb)
  2. Peak-based: TF motifs vs ATAC peaks → peak-to-gene links

Usage:
    python skeleton_build.py --dataset op
    python skeleton_build.py --dataset ibd_cd
"""

import argparse
import os
import sys
import anndata as ad
import networkx as nx
import pandas as pd
import scglue

env = os.environ

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", type=str, required=True, help="Dataset name (e.g. op, ibd_cd)")
args = parser.parse_args()

DATASET = args.dataset

par = {
    "rna":             f"{env['INFERENCE_DIR']}/{DATASET}_rna.h5ad",
    "atac":            f"{env['INFERENCE_DIR']}/{DATASET}_atac.h5ad",
    "annotation_file": f"{env['RESOURCES_DIR']}/supp_data/gencode.v47.annotation.gtf.gz",
    "motif_encode":    f"{env['RESOURCES_DIR']}/supp_data/databases/scglue/ENCODE-TF-ChIP-hg38.bed.gz",
    "motif_jaspar":    f"{env['RESOURCES_DIR']}/supp_data/databases/scglue/JASPAR2022-hg38.bed.gz",
    "temp_dir":        f"{env['RESULTS_DIR']}/experiment/skeleton/{DATASET}/tmp",
    "out_skeleton":    f"{env['RESULTS_DIR']}/experiment/skeleton/{DATASET}/skeleton.csv",
    "flank_length":    1000,
    "extend_range":    150_000,
}


def preprocess(par):
    """Annotate RNA var with genomic coords; build peak-to-gene guidance graph."""
    import scanpy as sc
    print("Reading RNA and ATAC ...", flush=True)
    rna  = ad.read_h5ad(par["rna"])
    atac = ad.read_h5ad(par["atac"])

    scglue.data.get_gene_annotation(rna, gtf=par["annotation_file"], gtf_by="gene_name")
    rna = rna[:, ~rna.var["chrom"].isna()]

    # scglue requires highly_variable column; compute if absent
    if "highly_variable" not in rna.var.columns:
        sc.pp.highly_variable_genes(rna, flavor="seurat_v3", n_top_genes=min(3000, rna.n_vars))
        rna.var["highly_variable"] = rna.var["highly_variable"].fillna(False)

    split = atac.var_names.str.split(r"[:-]")
    atac.var["chrom"]       = split.map(lambda x: x[0])
    atac.var["chromStart"]  = split.map(lambda x: x[1]).astype(int)
    atac.var["chromEnd"]    = split.map(lambda x: x[2]).astype(int)

    guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac, extend_range=par["extend_range"])
    scglue.graph.check_graph(guidance, [rna, atac])

    # Save annotated objects for downstream steps
    str_cols = [c for c in rna.var.columns if rna.var[c].dtype == object]
    rna.var[str_cols] = rna.var[str_cols].astype(str)

    rna.write(f"{par['temp_dir']}/rna.h5ad")
    atac.write(f"{par['temp_dir']}/atac.h5ad")
    nx.write_graphml(guidance, f"{par['temp_dir']}/guidance.graphml.gz")

    # Extract forward peak→gene edges
    fwd = guidance.edge_subgraph(
        e for e, attr in dict(guidance.edges).items() if attr["type"] == "fwd"
    ).copy()
    peak2gene = pd.DataFrame(
        [(e[1], e[0]) for e in fwd.edges], columns=["peak", "target"]
    )
    peak2gene.to_csv(f"{par['temp_dir']}/peak2gene.csv", index=False)
    print(f"peak2gene: {len(peak2gene)} entries", flush=True)


def build_skeleton_promoter(par, motif_file, label):
    """TF motifs vs TSS flanks → TF→gene edges."""
    rna = ad.read_h5ad(f"{par['temp_dir']}/rna.h5ad")
    flank_bed = scglue.genomics.Bed(rna.var).strand_specific_start_site().expand(
        par["flank_length"], par["flank_length"]
    )
    motif_bed = scglue.genomics.read_bed(motif_file)
    flank2tf  = scglue.genomics.window_graph(flank_bed, motif_bed, 0, right_sorted=True)
    df = pd.DataFrame(
        [(e[1], e[0]) for e in flank2tf.edges], columns=["source", "target"]
    )
    out = f"{par['temp_dir']}/skeleton_promoter_{label}.csv"
    df.to_csv(out, index=False)
    print(f"Promoter skeleton ({label}): {len(df)} edges", flush=True)
    return df


def build_skeleton_peak(par, motif_file, label):
    """TF motifs vs ATAC peaks, then bridged via peak→gene links."""
    atac      = ad.read_h5ad(f"{par['temp_dir']}/atac.h5ad")
    peak2gene = pd.read_csv(f"{par['temp_dir']}/peak2gene.csv")
    motif_bed = scglue.genomics.read_bed(motif_file)

    peak_bed  = scglue.genomics.Bed(atac.var)
    peak2tf   = scglue.genomics.window_graph(peak_bed, motif_bed, 0, right_sorted=True)
    df_peak2tf = pd.DataFrame(
        [(e[1], e[0]) for e in peak2tf.edges], columns=["source", "peak"]
    )
    # Merge: TF→peak→gene  →  TF→gene
    tf2gene = peak2gene.merge(df_peak2tf, on="peak", how="inner")[["source", "target"]].drop_duplicates()
    out = f"{par['temp_dir']}/skeleton_peak_{label}.csv"
    tf2gene.to_csv(out, index=False)
    print(f"Peak skeleton ({label}): {len(tf2gene)} edges", flush=True)
    return tf2gene


if __name__ == "__main__":
    os.makedirs(par["temp_dir"], exist_ok=True)
    os.makedirs(os.path.dirname(par["out_skeleton"]), exist_ok=True)

    # Step 1: preprocess (annotate genes, build peak→gene map)
    preprocess(par)

    motif_files = {"encode": par["motif_encode"], "jaspar": par["motif_jaspar"]}

    # Step 2: promoter-based skeletons
    promoter_parts = [build_skeleton_promoter(par, f, label) for label, f in motif_files.items()]
    skeleton_promoter = pd.concat(promoter_parts).drop_duplicates()

    # Step 3: peak-based skeletons
    peak_parts = [build_skeleton_peak(par, f, label) for label, f in motif_files.items()]
    skeleton_peak = pd.concat(peak_parts).drop_duplicates()

    # Step 4: merge and save
    skeleton = pd.concat([skeleton_promoter, skeleton_peak]).drop_duplicates()
    skeleton = skeleton[["source", "target"]]
    print(f"Final skeleton: {len(skeleton)} edges, "
          f"{skeleton['source'].nunique()} TFs, {skeleton['target'].nunique()} genes", flush=True)
    skeleton.to_csv(par["out_skeleton"], index=False)
    print(f"Saved to {par['out_skeleton']}", flush=True)
