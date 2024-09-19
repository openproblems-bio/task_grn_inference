import os
import anndata
import numpy as np
import pandas as pd
import subprocess
import ast
import requests

# wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

# wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather


## VIASH START
par = {
  'multiomics_rna': 'resources/grn-benchmark/multiomics_rna_0.h5ad',
  "tf_all": 'resources/prior/tf_all.csv',
  'prediction': 'output/scenic_0_hvgs.csv',
  'temp_dir': 'output/scenic',
  'num_workers': 20,
  'max_n_links': 50000,
  'rank_threshold': 10000,
  'auc_threshold': 0.05,
  'nes_threshold': 3.0,
  'seed': "32",
  'normalize': False,
  'only_hvgs': True
}
## VIASH END

import sys
meta= {
  "resources_dir": 'src/utils/'
}
sys.path.append(meta["resources_dir"])
from util import process_data, process_links
par['normalize']=False
# Load scRNA-seq data
print('Reading data')
adata_rna = anndata.read_h5ad(par['multiomics_rna'])
process_data(adata_rna, par)

os.makedirs(par['temp_dir'], exist_ok=True)

databases = f"{par['temp_dir']}/databases/"
os.makedirs(databases, exist_ok=True)

par['motif_annotation'] = f'{databases}/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'
par['genes_vs_motifs_10k'] = f'{databases}/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
par['genes_vs_motifs_500'] = f'{databases}/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'

if not (os.path.exists(par['motif_annotation'])):
  print('downloading motif_annotation')
  response = requests.get("https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
  with open(par['motif_annotation'], "wb") as file:
      file.write(response.content)
if not (os.path.exists(par['genes_vs_motifs_10k'])):
  print('downloading genes_vs_motifs_10k')
  response = requests.get("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
  with open(par['genes_vs_motifs_10k'], "wb") as file:
      file.write(response.content)
if not (os.path.exists(par['genes_vs_motifs_500'])):
  print('downloading genes_vs_motifs_500')
  response = requests.get("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
  with open(par['genes_vs_motifs_500'], "wb") as file:
      file.write(response.content)
      
expr_mat_adjacencies =  os.path.join(par['temp_dir'], "expr_mat_adjacencies.tsv")
expression_data = os.path.join(par['temp_dir'], "expression_data.tsv")
regulons = f"{par['temp_dir']}/regulons.csv"


def run_grn(par):
  print('Run grn')
  command = [
      "pyscenic", "grn",
      "--num_workers", str(par['num_workers']),
      "--seed", str(par['seed']),
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
      "--num_workers", str(par['num_workers']),
      "--rank_threshold", str(par['rank_threshold']),
      "--auc_threshold", str(par['auc_threshold']),
      "--nes_threshold", str(par['nes_threshold']),
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


