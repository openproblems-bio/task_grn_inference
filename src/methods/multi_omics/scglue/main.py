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

def preprocess(rna, atac, par):
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
    
    guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
    
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
        rna, "NB", use_highly_variable=False,
        use_layer="counts", use_rep="X_pca"
    )
    scglue.models.configure_dataset(
        atac, "NB", use_highly_variable=False,
        use_rep="X_lsi"
    )
    

    glue = scglue.models.fit_SCGLUE(
        {"rna": rna, "atac": atac}, guidance,
        fit_kws={"directory": f"{par['temp_dir']}/glue"}
    )

    glue.save(f"{par['temp_dir']}/glue.dill")

    dx = scglue.models.integration_consistency(
        glue, {"rna": rna, "atac": atac}, guidance
    )

    rna.obsm["X_glue"] = glue.encode_data("rna", rna)
    atac.obsm["X_glue"] = glue.encode_data("atac", atac)
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
    

def run_grn(par):
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
        if attr["qval"] < 0.05
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

    #Construct the command 
    # command = ['pyscenic', 'grn', f"{par['temp_dir']}/rna.loom", 
    #            f"{par['temp_dir']}/tfs.txt", '-o', f"{par['temp_dir']}/draft_grn.csv", 
    #            '--seed', '0', '--num_workers', f"{par['num_workers']}", 
    #            '--cell_id_attribute', 'obs_id', '--gene_attribute', 'name']
    # print('Run grn')
    # result = subprocess.run(command,  check=True)

    # print("Output:")
    # print(result.stdout)
    # print("Error:")
    # print(result.stderr)

    # if result.returncode == 0:
    #     print("Command executed successfully")
    # else:
    #     print("Command failed with return code", result.returncode)


    print("Generate TF cis-regulatory ranking bridged by ATAC peaks", flush=True)
    peak_bed = scglue.genomics.Bed(atac.var.loc[peaks])
    peak2tf = scglue.genomics.window_graph(peak_bed, motif_bed, 0, right_sorted=True)
    peak2tf = peak2tf.edge_subgraph(e for e in peak2tf.edges if e[1] in tfs)

    gene2tf_rank_glue = scglue.genomics.cis_regulatory_ranking(
        gene2peak, peak2tf, genes, peaks, tfs,
        region_lens=atac.var.loc[peaks, "chromEnd"] - atac.var.loc[peaks, "chromStart"],
        random_state=0)

    flank_bed = scglue.genomics.Bed(rna.var.loc[genes]).strand_specific_start_site().expand(10000, 10000)
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
def prune_grn(par):
    # Construct the command 
    #TODO: be sure that obs_id is in obs and name is in var
    print("Run pscenic ctx", flush=True)
    # command = [
    #     "pyscenic", "ctx",
    #     f"{par['temp_dir']}/draft_grn.csv",
    #     f"{par['temp_dir']}/glue.genes_vs_tracks.rankings.feather",
    #     f"{par['temp_dir']}/supp.genes_vs_tracks.rankings.feather",
    #     "--annotations_fname", f"{par['temp_dir']}/ctx_annotation.tsv",
    #     "--expression_mtx_fname", f"{par['temp_dir']}/rna.loom",
    #     "--output", f"{par['temp_dir']}/pruned_grn.csv",
    #     "--rank_threshold", "10000",
    #     "--auc_threshold", "0.1",
    #     "--nes_threshold", "2",
    #     "--mask_dropouts",
    #     "--min_genes", "1",
    #     "--num_workers", f"{par['num_workers']}",
    #     "--cell_id_attribute", "obs_id",
    #     "--gene_attribute", "name"
    # ]
    
    par['genes_vs_motifs_500'] = 'output/scenic/databases/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
    par['genes_vs_motifs_10k'] =  'output/scenic/databases/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
    par['motif_annotation'] = 'output/scenic/databases/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'

    command = [
        "pyscenic", "ctx", 
        f"{par['temp_dir']}/draft_grn.csv",
        f"{par['temp_dir']}/glue.genes_vs_tracks.rankings.feather",
        f"{par['temp_dir']}/supp.genes_vs_tracks.rankings.feather",
        "--annotations_fname", f"{par['temp_dir']}/ctx_annotation.tsv",
        "--expression_mtx_fname", f"{par['temp_dir']}/rna.loom",
        "--output", f"{par['temp_dir']}/pruned_grn.csv",
        # "--top_n_targets", "100",
        # "--rank_threshold", "1500",
        # "--auc_threshold", "0.1",
        # "--nes_threshold", "0",
        "--mask_dropouts",
        "--min_genes", "1",
        "--num_workers", f"{par['num_workers']}",
        "--cell_id_attribute", "obs_id",
        "--gene_attribute", "name"
    ]

    result = subprocess.run(command, check=True)

    print("Output:")
    print(result.stdout)
    print("Error:")
    print(result.stderr)

    if result.returncode == 0:
        print("pyscenic ctx executed successfully")
    else:
        print("pyscenic ctx failed with return code", result.returncode)

# def download_prior(par):
#     # get gene annotation
#     response = requests.get("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz")
#     par['annotation_file'] = f"{par['temp_dir']}/gencode.v45.annotation.gtf.gz"
#     if response.status_code == 200:
#         # Write the content to a file
#         with open(par['annotation_file'], 'wb') as file:
#             file.write(response.content)
#         print(f"File downloaded and saved as gencode.v45.annotation.gtf.gz")
#     else:
#         print(f"Failed to download the gencode.v45.annotation.gtf.gz. Status code: {response.status_code}")


    
#     annotation_file


def main(par):
    import torch
    print("Is CUDA available:", torch.cuda.is_available())
    print("Number of GPUs:", torch.cuda.device_count())
    
    os.makedirs(par['temp_dir'], exist_ok=True)
    print('Reading input files', flush=True)
    rna = ad.read_h5ad(par['multiomics_rna'])
    atac = ad.read_h5ad(par['multiomics_atac'])

    # print('Preprocess data', flush=True)
    # preprocess(rna, atac, par)
    # print('Train a model', flush=True)
    # training(par)
    run_grn(par)
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

    return scglue_grn
    