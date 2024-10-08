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
    feature_embeddings = glue.encode_graph(guidance)
    feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
    rna.varm["X_glue"] = feature_embeddings.reindex(rna.var_names).to_numpy()
    atac.varm["X_glue"] = feature_embeddings.reindex(atac.var_names).to_numpy()
    
    rna.write(f"{par['rna-emb']}", compression="gzip")
    atac.write(f"{par['atac-emb']}", compression="gzip")
    nx.write_graphml(guidance, f"{par['guidance.graphml']}")
def peak_tf_gene_connections(par):
    ''' Infers gene2peak connections
    '''
    print('reload the data')
    rna = ad.read_h5ad(f"{par['temp_dir']}/rna-emb.h5ad")
    atac = ad.read_h5ad(f"{par['temp_dir']}/atac-emb.h5ad")
    guidance = nx.read_graphml(f"{par['temp_dir']}/guidance.graphml.gz")

    rna.var["name"] = rna.var_names
    atac.var["name"] = atac.var_names

    genes = rna.var.index
    peaks = atac.var.index
    features = pd.Index(np.concatenate([rna.var_names, atac.var_names]))
    feature_embeddings = np.concatenate([rna.varm["X_glue"], atac.varm["X_glue"]])
    print('Get the skeleton')

    skeleton = guidance.edge_subgraph(
        e for e, attr in dict(guidance.edges).items()
        if attr["type"] == "fwd"
    ).copy()
    print('reginf')
    reginf = scglue.genomics.regulatory_inference(
        features, feature_embeddings,
        skeleton=skeleton, random_state=0
    )
    print('gene2peak')
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
    print('this is the motif file: ', par['motif_file'])
    motif_bed = scglue.genomics.read_bed(par['motif_file']) 
    # motif_bed = motif_bed.iloc[:100000, :] #TODO: remove this
    # tfs = pd.Index(motif_bed["name"]).intersection(rna.var_names)

    print("Generate TF cis-regulatory ranking bridged by ATAC peaks", flush=True)
    peak_bed = scglue.genomics.Bed(atac.var.loc[peaks])
    peak2tf = scglue.genomics.window_graph(peak_bed, motif_bed, 0, right_sorted=True)
    # peak2tf = peak2tf.edge_subgraph(e for e in peak2tf.edges if e[1] in tfs)

    flank_bed = scglue.genomics.Bed(rna.var.loc[genes]).strand_specific_start_site().expand(500, 500)
    flank2tf = scglue.genomics.window_graph(flank_bed, motif_bed, 0, right_sorted=True)
        
    sources = []
    targets = []
    for e, attr in dict(gene2peak.edges).items():
        sources.append(e[0])
        targets.append(e[1])
    df = pd.DataFrame({'source': sources, 'target':targets})
    df.to_csv(par['gene2peak'])  

    sources = []
    targets = []
    for e, attr in dict(peak2tf.edges).items():
        sources.append(e[0])
        targets.append(e[1])
    df = pd.DataFrame({'source': sources, 'target':targets})
    df.to_csv(par['peak2tf'])   

    sources = []
    targets = []
    for e, attr in dict(flank2tf.edges).items():
        sources.append(e[0])
        targets.append(e[1])
    df = pd.DataFrame({'source': sources, 'target':targets})
    df.to_csv(par['flank2tf']) 

def merge_connections(par):
    
    gene2peak = pd.read_csv(par['gene2peak'], index_col=0)  
    gene2peak.columns = ['target', 'peak']

    peak2tf= pd.read_csv(par['peak2tf'], index_col=0)   
    peak2tf.columns = ['peak', 'source']

    flank2tf= pd.read_csv(par['flank2tf'], index_col=0) 
    flank2tf.columns = ['target', 'source']
    # merge gene2peak and peak2tf
    tf2gene = gene2peak.merge(peak2tf, on='peak', how='inner')[['source','target']].drop_duplicates()
    # merge flank2tf and tf2gene
    tf2gene = pd.concat([tf2gene, flank2tf], axis=0).drop_duplicates()

    tf2gene.to_csv(f"{par['tf2gene']}")  

if __name__ == '__main__':
    par = {
        'multiomics_atac': f"resources/grn-benchmark/multiomics_atac.h5ad",
        'multiomics_rna': f"resources/grn-benchmark/multiomics_rna.h5ad",
        'annotation_file': f"output/db/gencode.v45.annotation.gtf.gz",
        # 'motif_file':   'output/db/ENCODE-TF-ChIP-hg38.bed.gz',
        'motif_file':   'output/db/jaspar_encode.bed.gz',
        'temp_dir': 'output/skeleton',
        'extend_range': 150000,
        'tf2gene': 'output/skeleton/tf2gene.csv'
    }
    print(par)
    os.makedirs(par['temp_dir'], exist_ok=True)
    par['rna-emb'] = f"{par['temp_dir']}/rna-emb.h5ad"
    par['atac-emb'] = f"{par['temp_dir']}/atac-emb.h5ad"
    par['guidance.graphml'] = f"{par['temp_dir']}/guidance.graphml.gz"
    
    par['gene2peak'] = f"{par['temp_dir']}/gene2peak.csv"
    par['peak2tf'] = f"{par['temp_dir']}/peak2tf.csv"
    par['flank2tf'] = f"{par['temp_dir']}/flank2tf.csv"

    # ---- simplify 
    if False:
        multiomics_atac = ad.read_h5ad(par['multiomics_atac'])
        multiomics_atac = multiomics_atac[:, :10000]

        par['multiomics_atac'] = f"{par['temp_dir']}/multiomics_atac.h5ad"
        multiomics_atac.write(par['multiomics_atac'])

    # ----- actual runs
    # print('------- preprocess ---------')
    # preprocess(par)
    # print('------- training ---------')
    # training(par)
    print('------- peak_tf_gene_connections ---------')
    peak_tf_gene_connections(par)
    print('------- merge_connections ---------')
    merge_connections(par)
    
    




