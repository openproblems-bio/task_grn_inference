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

## VIASH START
par = {
        'atac': f"resources/grn_benchmark/inference_data/op_atac.h5ad",
        'rna': f"resources/grn_benchmark/inference_data/op_rna.h5ad",
        'annotation_file': f"resources/supp_data/gencode.v45.annotation.gtf.gz",
        'temp_dir': 'output/skeleton',
        'extend_range': 150000,
        'flank_length': 1000,
        'skeleton': 'resources/grn_benchmark/prior/skeleton.csv',
        'motif_dataset_encode': 'resources/supp_data/databases/scglue/ENCODE-TF-ChIP-hg38.bed.gz',
        'motif_dataset_jaspar': 'resources/supp_data/databases/scglue/JASPAR2022-hg38.bed.gz'
    }
## VIASH END
def preprocess(par):
    print('Reading input files', flush=True)
    rna = ad.read_h5ad(par['rna'])
    atac = ad.read_h5ad(par['atac'])

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
    ]
    rna.var[column_names] = rna.var[column_names].astype(str)

    rna.write(f"{par['temp_dir']}/rna.h5ad")
    atac.write(f"{par['temp_dir']}/atac.h5ad")
    nx.write_graphml(guidance, f"{par['temp_dir']}/guidance.graphml.gz")

    guidance = guidance.edge_subgraph(
        e for e, attr in dict(guidance.edges).items()
        if attr["type"] == "fwd"
    ).copy()
    sources = []
    targets = []
    for e, attr in dict(guidance.edges).items():
        sources.append(e[1])
        targets.append(e[0])
    df = pd.DataFrame({'source': sources, 'target':targets})
    df.to_csv(par['peak2gene.csv'])

def get_flank_bed(par):
    rna = ad.read_h5ad(par['rna'])

    scglue.data.get_gene_annotation(
        rna, gtf=par['annotation_file'],
        gtf_by="gene_name"
    )
    rna = rna[:, ~rna.var.chrom.isna()]

    flank_bed = scglue.genomics.Bed(rna.var).strand_specific_start_site().expand(par['flank_length'], par['flank_length'])
    print(flank_bed.head(5))
    return flank_bed

def skeleton_promotor(par):
    '''Creates promotor based skeleton using TF motif data and TSS flank'''
    flank_bed = get_flank_bed(par)
    
    motif_bed = scglue.genomics.read_bed(par['motif_file']) 
    
    flank2tf = scglue.genomics.window_graph(flank_bed, motif_bed, 0, right_sorted=True)
    sources = []
    targets = []
    for e, attr in dict(flank2tf.edges).items():
        sources.append(e[1])
        targets.append(e[0])
    df = pd.DataFrame({'source': sources, 'target':targets})
    df.to_csv(par['skeleton_promotor_file'])

def skeleton_peak(par):
    '''Creates peak based skeleton using TF motif data'''
    atac = ad.read_h5ad(f"{par['temp_dir']}/atac.h5ad")
    
    print('this is the motif file: ', par['motif_file'])
    motif_bed = scglue.genomics.read_bed(par['motif_file']) 

    print("Generate TF cis-regulatory ranking bridged by ATAC peaks", flush=True)
    peak_bed = scglue.genomics.Bed(atac.var)
    peak2tf = scglue.genomics.window_graph(peak_bed, motif_bed, 0, right_sorted=True)
    peak2tf = peak2tf.edge_subgraph(e for e in peak2tf.edges if e[1] in tfs)

    sources = []
    targets = []
    for e, attr in dict(peak2tf.edges).items():
        sources.append(e[1])
        targets.append(e[0])
    peak2tf = pd.DataFrame({'source': sources, 'target':targets})

    # merge peak2tf with peak2gene
    peak2tf.columns = ['peak', 'source']

    peak2gene = pd.read_csv('output/skeleton/peak2gene.csv', index_col=0)  
    peak2gene.columns = ['peak', 'target']

    tf2gene = peak2gene.merge(peak2tf, on='peak', how='inner')[['source','target']].drop_duplicates() #TODO: fix this by adding peak


    tf2gene.to_csv(par['skeleton_peak_file']) 

if __name__ == '__main__':
    
    print(par)
    os.makedirs(par['temp_dir'], exist_ok=True)
    par['rna-emb'] = f"{par['temp_dir']}/rna-emb.h5ad"
    par['atac-emb'] = f"{par['temp_dir']}/atac-emb.h5ad"
    par['guidance.graphml'] = f"{par['temp_dir']}/guidance.graphml.gz"
    par['peak2gene'] = f"{par['temp_dir']}/peak2gene.csv"
    par['flank2gene'] = f"{par['temp_dir']}/flank2gene.csv"


    # ----- connect rna to atac -> peak2gene connections 
    if True:
        print('------- preprocess ---------')
        preprocess(par)


    print('------- promotor based skeleton for different motif files ---------')
    names = ['encode', 'jaspar']
    motif_files = [par['motif_dataset_jaspar'], par['motif_dataset_encode']]
    if True:
        for i, motif_file in enumerate(motif_files):
            par['skeleton_promotor_file'] = f"{par['temp_dir']}/skeleton_{names[i]}_promotor.csv"
            par['motif_file'] = motif_file
            skeleton_promotor(par)
        # - merge them 
        for i in range(len(names)):
            df = pd.read_csv(f"{par['temp_dir']}/skeleton_{names[i]}_promotor.csv", index_col=0)
            print(df.shape)
            if i ==0 :
                skeleton = df 
            else:
                skeleton = pd.concat([df, skeleton], axis=0).drop_duplicates()
            print(skeleton.shape)
        skeleton.to_csv(f"{par['temp_dir']}/skeleton_promotor.csv")
    print('------- peak based skeleton for different motif files ---------')
    if True:
        for i, motif_file in enumerate(motif_files):
            par['skeleton_peak_file'] = f"{par['temp_dir']}/skeleton_peak_{names[i]}.csv"
            par['motif_file'] = motif_file
            skeleton_peak(par)
        # - merge them 
        print('merging peak2tf from different motifs')
        for i, name in enumerate(names):
            df = pd.read_csv(f"{par['temp_dir']}/skeleton_peak_{names[i]}.csv")
            print(df.source.nunique())
            print(df.target.nunique())
            if i ==0 :
                skeleton_peak = df 
            else:
                skeleton_peak = pd.concat([df, skeleton_peak], axis=0).drop_duplicates()
        skeleton_peak.to_csv(f"{par['temp_dir']}/skeleton_peak.csv")

    if True:
        print('------- mege peak based skeleton with promotor base skeleton ---------')
        # - read peak based and promotor based skeletons
        skeleton_peak = pd.read_csv(f"{par['temp_dir']}/skeleton_peak.csv")[['source', 'target']].drop_duplicates()
        print(len(skeleton_peak), skeleton_peak.source.nunique(), skeleton_peak.target.nunique())
        skeleton_promotor = pd.read_csv(f"{par['temp_dir']}/skeleton_promotor.csv")[['source', 'target']].drop_duplicates()
        print(len(skeleton_promotor), skeleton_promotor.source.nunique(), skeleton_promotor.target.nunique())

        # - merge and save 
        skeleton = pd.concat([skeleton_promotor, skeleton_peak], axis=0).drop_duplicates()
        print(len(skeleton), skeleton.source.nunique(), skeleton.target.nunique())
        skeleton.to_csv(par['skeleton'])
    




