
import sys 
import os
## VIASH START
par = {
  'multiomics_rna': 'resources_test/grn-benchmark/multiomics_rna.h5ad',
  'multiomics_atac': 'resources_test/grn-benchmark/multiomics_atac.h5ad',
  'temp_dir': 'output/scenicplus',
  'prediction': 'output/prediction.csv',
  'qc': False,
  'num_workers': 4,
  'scplus_mdata': 'output/scenicplus/scplus_mdata.h5mu',
  'cell_topic': 'output/scenicplus/cell_topic.csv',
  'grn_extended': 'output/scenicplus/grn_extended.csv'
}
## VIASH END

import argparse
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--multiomics_rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--multiomics_atac', type=str, help='Path to the multiomics atac file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
args = parser.parse_args()

if args.multiomics_rna:
    par['multiomics_rna'] = args.multiomics_rna
if args.multiomics_atac:
    par['multiomics_atac'] = args.multiomics_atac
if args.prediction:
    par['prediction'] = args.prediction
if args.tf_all:
    par['tf_all'] = args.tf_all
if args.num_workers:
    par['num_workers'] = args.num_workers


if args.resources_dir:
    meta = {}
    meta['resources_dir'] = args.resources_dir  
par['num_workers'] = int(par['num_workers'])
print(par)

sys.path.append(meta["resources_dir"])
from main import * 


def print_memory_usage():
    import psutil

    process = psutil.Process(os.getpid())
    mem_info = process.memory_info().rss / (1024 * 1024)  # Convert to MB
    print(f"Memory usage: {mem_info:.2f} MB")


def main(par):

    par['cistopic_object'] = f'{par["temp_dir"]}/cistopic_object.pkl'
    par['DB_PATH'] = os.path.join('output', 'db')
    os.makedirs(par['DB_PATH'], exist_ok=True)
    par['motif_annotation'] = os.path.join(par['DB_PATH'], 'motifs-v10-nr.hgnc-m0.00001-o0.0.tbl')
    par['RANKINGS_DB_PATH'] = os.path.join(par['DB_PATH'], 'db.regions_vs_motifs.rankings.feather')
    par['SCORES_DB_PATH']= os.path.join(par['DB_PATH'], 'db.regions_vs_motifs.scores.feather')
    par['blacklist_path'] = os.path.join(par['DB_PATH'], 'hg38-blacklist.v2.bed')
    par['atac_dir'] = os.path.join(par['temp_dir'], 'atac')
    par['fragments_dict'] = os.path.join(par['temp_dir'], 'fragments_dict.json')
    par['MALLET_PATH'] = os.path.join(par['temp_dir'], 'Mallet-202108', 'bin', 'mallet')
    os.makedirs(par['atac_dir'], exist_ok=True)

    # print('------- download_databases -------')
    # download_databases(par)
    # print_memory_usage()
    # print('------- process_peak -------')
    # process_peak(par)
    # print_memory_usage()
    # print('------- run_cistopic -------')
    # run_cistopic(par)
    # print_memory_usage()
    # print('------- process_topics -------')
    # process_topics(par)
    # print_memory_usage()
    # print('------- preprocess_rna -------')
    # preprocess_rna(par)
    # print_memory_usage()
    # print('------- snakemake_pipeline -------')
    snakemake_pipeline(par)
    print_memory_usage()
    print('------- post_process -------')
    post_process(par)
if __name__ == '__main__':
    main(par)
