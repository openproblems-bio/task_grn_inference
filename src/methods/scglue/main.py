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
    

def cis_inference(par):
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

    # Construct the command 
    command = (
        f"pyscenic grn {par['temp_dir']}/rna.loom {par['temp_dir']}/tfs.txt "
        f"-o {par['temp_dir']}/draft_grn.csv --seed 0 --num_workers {par['num_workers']} "
        "--cell_id_attribute obs_id --gene_attribute name"
    )

    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    print("Output:")
    print(result.stdout)
    print("Error:")
    print(result.stderr)

    if result.returncode == 0:
        print("Command executed successfully")
    else:
        print("Command failed with return code", result.returncode)


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

    ### Prune coexpression network using cis-regulatory ranking

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

    # Construct the command 
    #TODO: be sure that obs_id is in obs and name is in var
    print("Run pscenic ctx", flush=True)
    command = (
        f" pyscenic ctx {par['temp_dir']}/draft_grn.csv {par['temp_dir']}/glue.genes_vs_tracks.rankings.feather "
        f" {par['temp_dir']}/supp.genes_vs_tracks.rankings.feather  --annotations_fname {par['temp_dir']}/ctx_annotation.tsv "
        f" --expression_mtx_fname {par['temp_dir']}/rna.loom --output {par['temp_dir']}/pruned_grn.csv "
        f" --rank_threshold 500 --min_genes 1  --num_workers {par['num_workers']} "
        " --cell_id_attribute obs_id --gene_attribute name"
    )

    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    print("Output:")
    print(result.stdout)
    print("Error:")
    print(result.stderr)

    if result.returncode == 0:
        print("pyscenic ctx executed successfully")
    else:
        print("pyscenic ctx failed with return code", result.returncode)

def main(par):
    
    
    os.makedirs(par['temp_dir'], exist_ok=True)
    print('Reading input files', flush=True)
    rna = ad.read_h5ad(par['multiomics_rna'])
    atac = ad.read_h5ad(par['multiomics_atac'])

    print('Preprocess data', flush=True)
    preprocess(rna, atac, par)
    print('Train a model', flush=True)
    training(par)
    cis_inference(par)
    print('Curate predictions', flush=True)
    df = pd.read_csv(
        f"{par['temp_dir']}/pruned_grn.csv", header=None, skiprows=3,
        usecols=[0, 8], names=["TF", "targets"]
    )

    tfs_list = []
    target_list = []
    weight_list = []
    for i, (tf, targets) in df.iterrows():
        for target, weight in literal_eval(targets):
            tfs_list.append(tf)
            target_list.append(target)
            weight_list.append(weight)
    scglue_grn = pd.DataFrame(np.stack([tfs_list, target_list, weight_list], axis=1), columns=['source','target','weight'])
    scglue_grn.weight = scglue_grn.weight.astype(float)
    scglue_grn = scglue_grn.drop_duplicates().reset_index(drop=True)

    
    # scglue_grn = pd.DataFrame(
    # data = {'source':['tf1'], 'target':['g1'], 'weight':[1]}
    # )

    return scglue_grn
    