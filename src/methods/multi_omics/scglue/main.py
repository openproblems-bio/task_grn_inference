import sys
import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
import os 
import subprocess
import pandas as pd
import numpy as np
from ast import literal_eval
import requests
import torch
def download_annotation(par):
    
    if not os.path.exists(par['annotation_file']):
        print("Downloading prior started")
    
        response = requests.get("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz")
        
        if response.status_code == 200:
            with open(par['annotation_file'], 'wb') as file:
                file.write(response.content)
            print(f"File downloaded and saved to {par['annotation_file']}")
        else:
            print(f"Failed to download the gencode.v45.annotation.gtf.gz. Status code: {response.status_code}")
        print("Downloading prior ended")
def download_motifs(par):
    # get gene annotation
    
    if not os.path.exists(par['motif_file']):
        tag = par['motif_file'].split('/')[-1]
        print("Downloading motif started")
        response = requests.get(f"http://download.gao-lab.org/GLUE/cisreg/{tag}")
        
        if response.status_code == 200:
            with open(par['motif_file'], 'wb') as file:
                file.write(response.content)
            print(f"File downloaded and saved to {par['motif_file']}")
        else:
            print(f"Failed to download the motif file. Status code: {response.status_code}")
        print("Downloading motif ended")
def preprocess(par):
    print('Reading input files', flush=True)
    rna = ad.read_h5ad(par['multiomics_rna'])
    atac = ad.read_h5ad(par['multiomics_atac'])
    
    rna.layers["counts"] = rna.X.copy()
    sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")
    sc.pp.normalize_total(rna)
    sc.pp.log1p(rna)
    sc.pp.scale(rna)
    sc.tl.pca(rna, n_comps=100, svd_solver="auto")
    sc.pp.neighbors(rna, metric="cosine")
    sc.tl.umap(rna)
    print('step 1 completed')
    
    scglue.data.lsi(atac, n_components=100, n_iter=15)
    sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
    sc.tl.umap(atac)
    print('step 2 completed')
    
    scglue.data.get_gene_annotation(
        rna, gtf=par['annotation_file'],
        gtf_by="gene_name"
    )
    
    rna = rna[:, ~rna.var.chrom.isna()]
    
    split = atac.var_names.str.split(r"[:-]")
    atac.var["chrom"] = split.map(lambda x: x[0])
    atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
    atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
    
    guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac, extend_range=par['extend_range'])
    
    scglue.graph.check_graph(guidance, [rna, atac])
    
    column_names = [
        "chrom",
        "gene_type",
        "gene_id",
        "hgnc_id",
        "havana_gene",
        "tag",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "itemRgb",
        "blockCount",
        "blockSizes",
        "blockStarts",
        "artif_dupl",
        "highly_variable_rank"
    ]
    rna.var[column_names] = rna.var[column_names].astype(str)

    rna.write(f"{par['temp_dir']}/rna.h5ad")
    atac.write(f"{par['temp_dir']}/atac.h5ad")
    nx.write_graphml(guidance, f"{par['temp_dir']}/guidance.graphml.gz")

def training(par):
    os.makedirs(f"{par['temp_dir']}/glue", exist_ok=True)
    rna = ad.read_h5ad(f"{par['temp_dir']}/rna.h5ad")
    atac = ad.read_h5ad(f"{par['temp_dir']}/atac.h5ad")
    guidance = nx.read_graphml(f"{par['temp_dir']}/guidance.graphml.gz")
    scglue.models.configure_dataset(
        rna, "NB", use_highly_variable=True,
        use_layer="counts", use_rep="X_pca", use_batch='donor_id', use_cell_type='cell_type'
    )
    scglue.models.configure_dataset(
        atac, "NB", use_highly_variable=True,
        use_rep="X_lsi", use_batch='donor_id', use_cell_type='cell_type'
    )
    if False:
        guidance_hvf = guidance.subgraph(chain(
            rna.var.query("highly_variable").index,
            atac.var.query("highly_variable").index
        )).copy()

    glue = scglue.models.fit_SCGLUE(
        {"rna": rna, "atac": atac}, guidance,
        fit_kws={"directory": f"{par['temp_dir']}/glue"}
    )

    glue.save(f"{par['temp_dir']}/glue.dill")

    if True: # consistency score
        dx = scglue.models.integration_consistency(
            glue, {"rna": rna, "atac": atac}, guidance
        )
        dx.to_csv(f"{par['temp_dir']}/consistency_scores.csv")

    rna.obsm["X_glue"] = glue.encode_data("rna", rna)
    atac.obsm["X_glue"] = glue.encode_data("atac", atac)
    if False:
        combined = ad.concat([rna, atac])
        sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
        sc.tl.umap(combined)
    feature_embeddings = glue.encode_graph(guidance)
    feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
    rna.varm["X_glue"] = feature_embeddings.reindex(rna.var_names).to_numpy()
    atac.varm["X_glue"] = feature_embeddings.reindex(atac.var_names).to_numpy()

    rna.write(f"{par['temp_dir']}/rna-emb.h5ad", compression="gzip")
    atac.write(f"{par['temp_dir']}/atac-emb.h5ad", compression="gzip")
    nx.write_graphml(guidance, f"{par['temp_dir']}/guidance.graphml.gz")
    
def create_prior(par):
    ''' Infers gene2peak connections
    '''
    rna = ad.read_h5ad(f"{par['temp_dir']}/rna-emb.h5ad")
    atac = ad.read_h5ad(f"{par['temp_dir']}/atac-emb.h5ad")
    guidance = nx.read_graphml(f"{par['temp_dir']}/guidance.graphml.gz")

    rna.var["name"] = rna.var_names
    atac.var["name"] = atac.var_names

    genes = rna.var.index
    peaks = atac.var.index

    features = pd.Index(np.concatenate([rna.var_names, atac.var_names]))
    feature_embeddings = np.concatenate([rna.varm["X_glue"], atac.varm["X_glue"]])

    skeleton = guidance.edge_subgraph(
        e for e, attr in dict(guidance.edges).items()
        if attr["type"] == "fwd"
    ).copy()

    reginf = scglue.genomics.regulatory_inference(
        features, feature_embeddings,
        skeleton=skeleton, random_state=0
    )

    gene2peak = reginf.edge_subgraph(
        e for e, attr in dict(reginf.edges).items()
        if attr["qval"] < 0.1
    )

    scglue.genomics.Bed(atac.var).write_bed(f"{par['temp_dir']}/peaks.bed", ncols=3)
    scglue.genomics.write_links(
        gene2peak,
        scglue.genomics.Bed(rna.var).strand_specific_start_site(),
        scglue.genomics.Bed(atac.var),
        f"{par['temp_dir']}/gene2peak.links", keep_attrs=["score"]
    )

    motif_bed = scglue.genomics.read_bed(par['motif_file']) ## http://download.gao-lab.org/GLUE/cisreg/JASPAR2022-hg38.bed.gz
    tfs = pd.Index(motif_bed["name"]).intersection(rna.var_names)
    rna[:, np.union1d(genes, tfs)].write_loom(f"{par['temp_dir']}/rna.loom")
    np.savetxt(f"{par['temp_dir']}/tfs.txt", tfs, fmt="%s")

    print("Generate TF cis-regulatory ranking bridged by ATAC peaks", flush=True)
    peak_bed = scglue.genomics.Bed(atac.var.loc[peaks])
    peak2tf = scglue.genomics.window_graph(peak_bed, motif_bed, 0, right_sorted=True)
    peak2tf = peak2tf.edge_subgraph(e for e in peak2tf.edges if e[1] in tfs)

    gene2tf_rank_glue = scglue.genomics.cis_regulatory_ranking(
        gene2peak, peak2tf, genes, peaks, tfs,
        region_lens=atac.var.loc[peaks, "chromEnd"] - atac.var.loc[peaks, "chromStart"],
        random_state=0)

    flank_bed = scglue.genomics.Bed(rna.var.loc[genes]).strand_specific_start_site().expand(500, 500)
    flank2tf = scglue.genomics.window_graph(flank_bed, motif_bed, 0, right_sorted=True)

    gene2flank = nx.Graph([(g, g) for g in genes])
    gene2tf_rank_supp = scglue.genomics.cis_regulatory_ranking(
        gene2flank, flank2tf, genes, genes, tfs,
        n_samples=0
    )

    ### Prepare data for pruning 
    print("Prepare data for pruning ")

    gene2tf_rank_glue.columns = gene2tf_rank_glue.columns + "_glue"
    gene2tf_rank_supp.columns = gene2tf_rank_supp.columns + "_supp"

    scglue.genomics.write_scenic_feather(gene2tf_rank_glue, f"{par['temp_dir']}/glue.genes_vs_tracks.rankings.feather")
    scglue.genomics.write_scenic_feather(gene2tf_rank_supp, f"{par['temp_dir']}/supp.genes_vs_tracks.rankings.feather")

    pd.concat([
        pd.DataFrame({
                "#motif_id": tfs + "_glue",
                "gene_name": tfs
            }),
        pd.DataFrame({
            "#motif_id": tfs + "_supp",
            "gene_name": tfs
        })
    ]).assign(
        motif_similarity_qvalue=0.0,
        orthologous_identity=1.0,
        description="placeholder"
    ).to_csv(f"{par['temp_dir']}/ctx_annotation.tsv", sep="\t", index=False)
def pyscenic_grn(par):
    command = ['pyscenic', 'grn', f"{par['temp_dir']}/rna.loom", 
            f"{par['temp_dir']}/tfs.txt", '-o', f"{par['temp_dir']}/draft_grn.csv", 
            '--seed', '0', '--num_workers', f"{par['num_workers']}", 
            '--cell_id_attribute', 'obs_id', '--gene_attribute', 'name']
    print('Run grn')
    result = subprocess.run(command,  check=True)

    print("Output:")
    print(result.stdout)
    print("Error:")
    print(result.stderr)

    if result.returncode == 0:
        print("Command executed successfully")
    else:
        print("Command failed with return code", result.returncode)
def prune_grn(par):
    # Construct the command 
    print(par)
    print("Run pyscenic ctx", flush=True)
    command = [
        "pyscenic", "ctx", 
        f"{par['temp_dir']}/draft_grn.csv",
        f"{par['temp_dir']}/glue.genes_vs_tracks.rankings.feather",
        f"{par['temp_dir']}/supp.genes_vs_tracks.rankings.feather",
        "--annotations_fname", f"{par['temp_dir']}/ctx_annotation.tsv",
        "--expression_mtx_fname", f"{par['temp_dir']}/rna.loom",
        "--output", f"{par['temp_dir']}/pruned_grn.csv",
        "--top_n_targets", str(par['top_n_targets']),
        "--rank_threshold", str(par['rank_threshold']),
        "--auc_threshold", "0",
        "--nes_threshold", "0", 
        "--min_genes", "1",
        "--thresholds", "0.5", "0.7", 
        "--top_n_regulators", "10", "50", "100",
        "--max_similarity_fdr", "0.1",
        "--num_workers", f"{par['num_workers']}",
        "--cell_id_attribute", "obs_id", # be sure that obs_id is in obs and name is in var
        "--gene_attribute", "name"
    ]
    
#       --thresholds, The first method to create the TF-modules based on the best targets for each transcription factor (default: 0.75 0.90).
#   --top_n_targets TOP_N_TARGETS [TOP_N_TARGETS ...]
#                         The second method is to select the top targets for a given TF. (default: 50)
# --top_n_regulators TOP_N_REGULATORS [TOP_N_REGULATORS ...]
#                         The alternative way to create the TF-modules is to select the best regulators for each gene. (default: 5 10 50)
# --max_similarity_fdr MAX_SIMILARITY_FDR
#                         Maximum FDR in motif similarity to use when annotating enriched motifs (default: 0.001)
# --auc_threshold AUC_THRESHOLD
#                         The threshold used for calculating the AUC of a feature as fraction of ranked genes (default: 0.05).
# --nes_threshold NES_THRESHOLD
#                         The Normalized Enrichment Score (NES) threshold for finding enriched features (default: 3.0).

    result = subprocess.run(command, check=True)
    print("Output:")
    print(result.stdout)
    print("Error:")
    print(result.stderr)

    if result.returncode == 0:
        print("pyscenic ctx executed successfully")
    else:
        print("pyscenic ctx failed with return code", result.returncode)


def main(par):
    
    print("Is CUDA available:", torch.cuda.is_available())
    print("Number of GPUs:", torch.cuda.device_count())
    
    from util import process_links
    # Load scRNA-seq data
    # os.makedirs(par['temp_dir'], exist_ok=True)
    # print('----- download_annotation ---- ', flush=True)
    # download_annotation(par)
    # print('----- download_motifs ---- ', flush=True)
    # download_motifs(par)
    # print('----- preprocess ---- ', flush=True)
    # preprocess(par)
    # print('----- training ---- ', flush=True)
    # training(par)
    # print('----- create_prior ---- ', flush=True)
    # create_prior(par)
    # print('----- pyscenic_grn ---- ', flush=True)
    # pyscenic_grn(par)
    # print('----- prune_grn ---- ', flush=True)
    prune_grn(par)
    print('Curate predictions', flush=True)
    pruned_grn = pd.read_csv(
        f"{par['temp_dir']}/pruned_grn.csv", header=None, skiprows=3,
        usecols=[0, 8], names=["TF", "targets"]
    )

    tfs_list = []
    target_list = []
    weight_list = []
    for i, (tf, targets) in pruned_grn.iterrows():
        for target, weight in literal_eval(targets):
            tfs_list.append(tf)
            target_list.append(target)
            weight_list.append(weight)
    scglue_grn = pd.DataFrame(np.stack([tfs_list, target_list, weight_list], axis=1), columns=['source','target','weight'])
    scglue_grn.weight = scglue_grn.weight.astype(float)
    scglue_grn = scglue_grn.drop_duplicates().reset_index(drop=True)
    
    # scglue_grn = process_links(scglue_grn, par)

    return scglue_grn
    