import os
import anndata
import numpy as np
import pandas as pd
import subprocess
import ast
import requests
import scipy.sparse as sp
import sys
import anndata as ad
import argparse


# wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

# wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather


## VIASH START
par = {
  'rna': 'resources/grn_benchmark/inference_data/replogle_rna.h5ad',
  "tf_all": 'resources/grn_benchmark/prior/tf_all.csv',
  'prediction': 'resources/results/replogle/replogle.scenic.scenic.prediction.h5ad',
  'temp_dir': 'output/scenic',
  'num_workers': 20,
  'max_n_links': 50000,
  'rank_threshold': 10000,
  'auc_threshold': 0.05,
  'nes_threshold': 3.0,
  'seed': "32",
  'normalize': False
}
## VIASH END

import argparse
argparser = argparse.ArgumentParser()
argparser.add_argument('--rna', type=str, help='Path to the input RNA data in h5ad format.')
argparser.add_argument('--prediction', type=str, help='Path to the output prediction in h5ad format.')
args = argparser.parse_args()
if args.rna is not None:
  par['rna'] = args.rna
if args.prediction is not None:
  par['prediction'] = args.prediction



try:
    sys.path.append(meta["resources_dir"])
    from util import process_links
    from helper import format_data, run_grn
except:
    meta = {
      'util_dir': 'src/utils',
      'helper_dir': 'src/methods/grnboost',
    }
    sys.path.append(meta["helper_dir"])
    from helper import format_data, run_grn
    sys.path.append(meta["util_dir"])
    from util import process_links
    


def download_prior(par):
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


def prune_grn(par):
  print('Run prune')
  
  command = [
      "pyscenic", "ctx",
      par['expr_mat_adjacencies'], par['genes_vs_motifs_500'], par['genes_vs_motifs_10k'],
      "--annotations_fname", par['motif_annotation'], 
      "--expression_mtx_fname", par['expression_data'],
      "--mode", "custom_multiprocessing",
      "--output", str(par['regulons']), 
      "--num_workers", str(par['num_workers']),
      "--rank_threshold", str(par['rank_threshold']),
      "--auc_threshold", str(par['auc_threshold']),
      "--nes_threshold", str(par['nes_threshold']),
  ]

  subprocess.run(command, check=True)
def format_grn(par):
    
  print('Format regulons')
  regulons_df = pd.read_csv(par['regulons'], index_col=0, skiprows=1)['TargetGenes'].reset_index().iloc[1:,:]

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
def main(par):
  databases = f"{par['temp_dir']}/databases/"
  os.makedirs(databases, exist_ok=True)
  os.makedirs(par['temp_dir'], exist_ok=True)

  par['motif_annotation'] = f'{databases}/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'
  par['genes_vs_motifs_10k'] = f'{databases}/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
  par['genes_vs_motifs_500'] = f'{databases}/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
  
  par['expr_mat_adjacencies'] =  os.path.join(par['temp_dir'], "expr_mat_adjacencies.tsv")
  par['expression_data'] = os.path.join(par['temp_dir'], "expression_data.tsv")
  par['regulons'] = f"{par['temp_dir']}/regulons.csv"
  
  download_prior(par)
  format_data(par)
  run_grn(par)
  prune_grn(par)
  network = format_grn(par)
  return network
  
if __name__=='__main__':
  dataset_id = ad.read_h5ad(par['rna'], backed='r').uns['dataset_id']
  net = main(par)
  os.makedirs(par['temp_dir'], exist_ok=True)

  net['weight'] = net['weight'].astype(str)
  output = ad.AnnData(X=None, uns={"method_id": "scenic", "dataset_id": dataset_id, "prediction": net[["source", "target", "weight"]]})
  output.write(par['prediction'])



