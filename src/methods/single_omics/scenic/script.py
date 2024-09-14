import os

import anndata
import numpy as np
import pandas as pd
# from contextlib_chdir import chdir 
import subprocess
# from urllib.request import urlretrieve
import zipfile
from arboreto.algo import grnboost2
from distributed import Client
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

out_dir = 'resources_local'
RANKINGS_DB_PATH = os.path.join(out_dir, 'cistarget-db', 'db.regions_vs_motifs.rankings.feather')
SCORES_DB_PATH = os.path.join(out_dir, 'cistarget-db', 'db.regions_vs_motifs.scores.feather')

# if not (os.path.exists(RANKINGS_DB_PATH) and os.path.exists(SCORES_DB_PATH)):

#     # Download create_cisTarget_databases
#     os.makedirs(os.path.join(out_dir, 'cistarget-db'), exist_ok=True)
#     if not os.path.exists(os.path.join(out_dir, 'cistarget-db', 'create_cisTarget_databases')):
#         with chdir(os.path.join(out_dir, 'cistarget-db')):
#             subprocess.run(['git', 'clone', 'https://github.com/aertslab/create_cisTarget_databases'])

#     # Download cluster-buster
#     if not os.path.exists(os.path.join(out_dir, 'cistarget-db', 'cbust')):
#         urlretrieve('https://resources.aertslab.org/cistarget/programs/cbust', os.path.join(out_dir, 'cistarget-db', 'cbust'))
#     subprocess.run(['chmod', 'a+x', os.path.join(out_dir, 'cistarget-db', 'cbust')])

#     # Download motif collection
#     if not os.path.exists(os.path.join(out_dir, 'cistarget-db', 'v10nr_clust_public')):
#         urlretrieve(
#             'https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip',
#             os.path.join(out_dir, 'cistarget-db', 'v10nr_clust_public.zip')
#         )
#         with zipfile.ZipFile(os.path.join(out_dir, 'cistarget-db', 'v10nr_clust_public.zip'), 'r') as zip_ref:
#             zip_ref.extractall(os.path.join(out_dir, 'cistarget-db'))

#     # Download chromosome sizes
#     if not os.path.exists(os.path.join(out_dir, 'cistarget-db', 'hg38.chrom.sizes')):
#         urlretrieve(
#             'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
#             os.path.join(out_dir, 'cistarget-db', 'hg38.chrom.sizes')
#         )

#     # Download reference genome
#     if not os.path.exists(os.path.join(out_dir, 'cistarget-db', 'hg38.fa')):
#         print('Downloading reference genome', flush=True)
#         if not os.path.exists(os.path.join(out_dir, 'cistarget-db', 'hg38.fa.gz')):
#             urlretrieve(
#                 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz',
#                 os.path.join(out_dir, 'cistarget-db', 'hg38.fa.gz')
#             )
#         with gzip.open(os.path.join(out_dir, 'cistarget-db', 'hg38.fa.gz'), 'rb') as f_in:
#             with open(os.path.join(out_dir, 'cistarget-db', 'hg38.fa'), 'wb') as f_out:
#                 shutil.copyfileobj(f_in, f_out)

#     # Prepare fasta from consensus regions
#     if not os.path.exists(os.path.join(out_dir, 'cistarget-db', 'hg38.with_1kb_bg_padding.fa')):
#         subprocess.run([
#             os.path.join(out_dir, 'cistarget-db', 'create_cisTarget_databases', 'create_fasta_with_padded_bg_from_bed.sh'),
#             os.path.join(out_dir, 'cistarget-db', 'hg38.fa'),
#             os.path.join(out_dir, 'cistarget-db', 'hg38.chrom.sizes'),
#             os.path.join(out_dir, 'consensus_peak_calling', 'consensus_regions.bed'),
#             os.path.join(out_dir, 'cistarget-db', 'hg38.with_1kb_bg_padding.fa'),
#             '1000',
#             'yes'
#         ])

#     # Create cistarget databases
#     with open(os.path.join(out_dir, 'cistarget-db', 'motifs.txt'), 'w') as f:
#         for filename in os.listdir(os.path.join(out_dir, 'cistarget-db', 'v10nr_clust_public', 'singletons')):
#             f.write(f'{filename}\n')
#     with chdir(os.path.join(out_dir, 'cistarget-db')):
#         subprocess.run([
#             'python',
#             os.path.join(out_dir, 'cistarget-db', 'create_cisTarget_databases', 'create_cistarget_motif_databases.py'),
#             '-f', os.path.join(out_dir, 'cistarget-db', 'hg38.with_1kb_bg_padding.fa'),
#             '-M', os.path.join(out_dir, 'cistarget-db', 'v10nr_clust_public', 'singletons'),
#             '-m', os.path.join(out_dir, 'cistarget-db', 'motifs.txt'),
#             '-c', os.path.join(out_dir, 'cistarget-db', 'cbust'),
#             '-o', 'db',
#             '--bgpadding', '1000',
#             '-t', str(par['num_workers'])
#         ])




# Load scRNA-seq data
adata_rna = anndata.read_h5ad(par['multiomics_rna'])
adata_rna = adata_rna[adata_rna.obs.donor_id=='donor_0',] #TODO: to go
gene_names = adata_rna.var.gene_ids.index.to_numpy()
X = adata_rna.X.toarray()

# Load list of putative TFs
# df = pd.read_csv(par["tf_all"], header=None, names=['gene_name'])
# tfs = set(list(df['gene_name']))
# tf_names = [gene_name for gene_name in gene_names if (gene_name in tfs)]

# format output
expression_data = os.path.join(par['temp_dir'], "expression_data.tsv")

pd.DataFrame(X, columns=gene_names).to_csv(expression_data, sep='\t', index=False)

expr_mat_adjacencies =  os.path.join(par['temp_dir'], "expr_mat_adjacencies.tsv")


print('Run grn')
command = [
    "pyscenic", "grn",
    "--num_workers", par['max_workers'],
    "--seed", par['seed'],
    "-o", expr_mat_adjacencies,
    expression_data,
    par['tf_all'],
    
]
subprocess.run(command, check=True)


print('Run prune')
regulons = f"{par['temp_dir']}/regulons.csv"
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

print('Format regulons')
regulons = pd.read_csv(regulons, index_col=0, skiprows=1)['TargetGenes'].reset_index().iloc[1:,:]

def format_it(df):
    values = df['TargetGenes'].values[0]
    genes = []
    weights = []
    for value in ast.literal_eval(values):
      genes.append(value[0])
      weights.append(value[1])
    return pd.DataFrame({'target':genes, 'weight':weights})
grn = regulons.groupby('index').apply(lambda df: format_it(df)).reset_index().rename(columns={'index':'source'})
network = grn[['source','target','weight']]


network.to_csv(par['prediction'], sep=',')

print('Finished.')


