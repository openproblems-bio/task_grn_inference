import os
import anndata
import numpy as np
import pandas as pd
import subprocess
import ast

# wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

# wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather


## VIASH START
par = {
  'multiomics_rna': 'resources/grn-benchmark/multiomics_rna.h5ad',
  "tf_all": 'resources/prior/tf_all.csv',
  'prediction': 'output/scenic/scenic.csv',
  'temp_dir': 'output/scenic',
  'max_workers': "10",
  'motif_annotation': 'resources_local/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl', 
  'genes_vs_motifs_10k': 'resources_local/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',
  'genes_vs_motifs_500': 'resources_local/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',
  'max_n_links': 50000,
  'seed': "32"
}
## VIASH END
os.makedirs(par['temp_dir'], exist_ok=True)

# Load list of putative TFs
# df = pd.read_csv(par["tf_all"], header=None, names=['gene_name'])
# tfs = set(list(df['gene_name']))
# tf_names = [gene_name for gene_name in gene_names if (gene_name in tfs)]


expr_mat_adjacencies =  os.path.join(par['temp_dir'], "expr_mat_adjacencies.tsv")
expression_data = os.path.join(par['temp_dir'], "expression_data.tsv")
regulons = f"{par['temp_dir']}/regulons.csv"

def format_data(par):
  print('Read data')
  adata_rna = anndata.read_h5ad(par['multiomics_rna'])  
  gene_names = adata_rna.var_names
  pd.DataFrame(adata_rna.X.todense(), columns=gene_names).to_csv(expression_data, sep='\t', index=False)
  

def run_grn(par):
  print('Run grn')
  command = [
      "pyscenic", "grn",
      "--num_workers", par['max_workers'],
      "--seed", par['seed'],
      "-o", expr_mat_adjacencies,
      "--method", "grnboost2", 
      expression_data,
      par['tf_all']
  ]
  subprocess.run(command, check=True)

def prune_grn(par):
  print('Run prune')
  
  command = [
      "pyscenic", "ctx",
      expr_mat_adjacencies, par['genes_vs_motifs_500'], par['genes_vs_motifs_10k'],
      "--annotations_fname", par['motif_annotation'], 
      "--expression_mtx_fname", expression_data,
      "--mode", "custom_multiprocessing",
      "--output", regulons, 
      "--num_workers", par['max_workers']
  ]

  subprocess.run(command, check=True)
def format_grn(par):
    
  print('Format regulons')
  regulons_df = pd.read_csv(regulons, index_col=0, skiprows=1)['TargetGenes'].reset_index().iloc[1:,:]

  def format_it(df):
      values = df['TargetGenes'].values[0]
      genes = []
      weights = []
      for value in ast.literal_eval(values):
        genes.append(value[0])
        weights.append(value[1])
      return pd.DataFrame({'target':genes, 'weight':weights})
  grn = regulons_df.groupby('index').apply(lambda df: format_it(df)).reset_index().rename(columns={'index':'source'})
  network = grn[['source','target','weight']]
  return network

format_data(par)
run_grn(par)
prune_grn(par)
network = format_grn(par)
network.to_csv(par['prediction'], sep=',')

print('Finished.')

