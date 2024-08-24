import os
import sys
import yaml
import pickle
import contextlib
import hashlib
import requests
import subprocess
from urllib.request import urlretrieve

import numpy as np
import scanpy as sc
import pandas as pd
import anndata
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
from scenicplus.wrappers.run_pycistarget import run_pycistarget


## VIASH START
par = {
  'multiomics_rna': 'resources/grn-benchmark/multiomics_rna.h5ad',
  'multiomics_atac': 'resources/grn-benchmark/multiomics_atac.h5ad',
  'temp_dir': 'output/scenicplus',
  'prediction': 'output/prediction.csv',
}
## VIASH END

work_dir = par['temp_dir']
par['cistopic_out'] = f'{work_dir}/cistopic_out'
par['cistopic_object'] = os.path.join(par['cistopic_out'], f'cistopic_object_with_model.pkl')
os.makedirs(os.path.join(work_dir, 'scRNA'), exist_ok=True)

# Download databases
DB_PATH = os.path.join(work_dir, 'db')
os.makedirs(DB_PATH, exist_ok=True)
def download(url: str, filepath: str) -> None:
    if os.path.exists(filepath):
        return
    print(f'Download {url}...')
    urlretrieve(url, filepath)
def download_and_checksum(url: str, filepath: str, digest: str) -> None:
    download(url, filepath)
    #with open(filepath, 'rb') as f:
    #    file_hash = hashlib.file_digest(f, 'sha1')
    #if file_hash.hexdigest() != digest:
    #    os.remove(filepath)
    #print(file_hash.hexdigest(), digest)
    #assert file_hash.hexdigest() == digest
def download_checksum(url: str, filepath: str) -> str:
    if not os.path.exists(filepath):
        response = requests.get(url)
        with open(filepath, 'w') as f:
            f.write(response.text)
    with open(filepath, 'r') as f:
        s = f.read()
    return s.split()[0]
url = 'https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather.sha1sum.txt'
digest = download_checksum(url, os.path.join(DB_PATH, 'hg38_screen_v10_clust.regions_vs_motifs.rankings.feather.sha1sum.txt'))
url = 'https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather'
download_and_checksum(url, os.path.join(DB_PATH, 'hg38_screen_v10_clust.regions_vs_motifs.rankings.feather'), digest)
url = 'https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.scores.feather.sha1sum.txt'
digest = download_checksum(url, os.path.join(DB_PATH, 'hg38_screen_v10_clust.regions_vs_motifs.scores.feather.sha1sum.txt'))
url = 'https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.scores.feather'
download_and_checksum(url, os.path.join(DB_PATH, 'hg38_screen_v10_clust.regions_vs_motifs.scores.feather'), digest)
url = 'https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/snapshots/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl'
download(url, os.path.join(DB_PATH, 'motifs-v10-nr.hgnc-m0.00001-o0.0.tbl'))

if not os.path.exists(os.path.join(work_dir, 'rna.h5ad')):

    # Load scRNA-seq data
    adata_rna = anndata.read_h5ad(par['multiomics_rna'])

    # Only keep cells with at least 200 expressed genes, and genes with at least 3 cells expressing them
    sc.pp.filter_cells(adata_rna, min_genes=200)
    sc.pp.filter_genes(adata_rna, min_cells=3)

    # Filter out doublets using scrublet
    sc.external.pp.scrublet(adata_rna)
    adata_rna = adata_rna[adata_rna.obs['predicted_doublet'] == False]

    # Filter based on mitochondrial and total counts
    adata_rna.var['mt'] = adata_rna.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata_rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Normalize data but keep track of original data
    adata_rna.raw = adata_rna
    sc.pp.normalize_total(adata_rna, target_sum=1e4)
    sc.pp.log1p(adata_rna)
    adata_rna.write_h5ad(os.path.join(work_dir, 'rna.h5ad'))

# Load candidate enhancer regions
with open(os.path.join(par['cistopic_out'], f'candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb') as f:
    region_bin_topics_otsu = pickle.load(f)
with open(os.path.join(par['cistopic_out'], f'candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb') as f:
    region_bin_topics_top3k = pickle.load(f)
with open(os.path.join(par['cistopic_out'], f'candidate_enhancers/markers_dict.pkl'), 'rb') as f:
    markers_dict = pickle.load(f)


# Convert to dictionary of pyrange objects
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}
region_sets['DARs'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')]
    region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

# Define rankings, score and motif annotation databases
rankings_db = os.path.join(DB_PATH, 'hg38_screen_v10_clust.regions_vs_motifs.rankings.feather')
scores_db =  os.path.join(DB_PATH, 'hg38_screen_v10_clust.regions_vs_motifs.scores.feather')
motif_annotation = os.path.join(DB_PATH, 'motifs-v10-nr.hgnc-m0.00001-o0.0.tbl')

# Init scenicplus pipeline
os.makedirs(os.path.join(work_dir, 'scplus_pipeline'), exist_ok=True)
os.makedirs(os.path.join(work_dir, 'scplus_pipeline', 'temp'), exist_ok=True)
subprocess.run(['scenicplus', 'init_snakemake', '--temp_dir', os.path.join(work_dir, 'scplus_pipeline')])

# Load pipeline settings
with open(os.path.join(work_dir, 'scplus_pipeline', 'Snakemake', 'config', 'config.yaml'), 'r') as f:
    settings = yaml.safe_load(f)

# Update settings: indicate locations of input files
settings['input_data']['cisTopic_obj_fname'] = par['cistopic_object']
settings['input_data']['GEX_anndata_fname'] = os.path.join(work_dir, 'rna.h5ad')
settings['input_data']['region_set_folder'] = os.path.join(par['cistopic_out'], 'region_sets')
settings['input_data']['ctx_db_fname'] = rankings_db
settings['input_data']['dem_db_fname'] = scores_db
settings['input_data']['path_to_motif_annotations'] = motif_annotation
settings['params_general']['temp_dir'] = os.path.join(work_dir, 'scplus_pipeline', 'temp')
settings['params_general']['n_cpu'] = 1

# Save pipeline settings
with open(os.path.join(work_dir, 'scplus_pipeline', 'Snakemake', 'config', 'config.yaml'), 'w') as f:
    yaml.dump(settings, f)

# TODO: from this line onward, the code is untested (could not run it locally due to excessive memory requirements)

# Run pipeline
with contextlib.chdir(os.path.join(work_dir, 'scplus_pipeline', 'Snakemake')):
    subprocess.run([
        'snakemake',
        '--cores', '1',
        #'--unlock'
    ])

# Make sure the file is properly formatted, and re-format it if needed
filepath = os.path.join(work_dir, 'tf_to_gene_adj.tsv')
shutil.copyfile(filepath, par['prediction'])


# cistopic_obj = pickle.load(os.path.join(par['cistopic_out'], f'cistopic_object_with_model.pkl'))
# # get cell topic association 
# cell_topic = cistopic_obj.selected_model.cell_topic.T
# cell_names = cistopic_obj.cell_data.obs_id.values
# cell_topic.index = cell_names
