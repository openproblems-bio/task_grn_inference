
import sys 
import os
import anndata as ad
import pandas as pd
import requests
## VIASH START
par = {
  'rna': 'resources_test/grn_benchmark/inference_data/op_rna.h5ad',
  'atac': 'resources_test/grn_benchmark/inference_data/op_atac.h5ad',
  'temp_dir': 'output/sp_new',
  'prediction': 'output/sp_new/prediction.h5ad',
  'qc': False,
  'num_workers': 20,
  'scplus_mdata': 'output/sp_new/scplus_mdata.h5mu',
  'cell_topic': 'output/sp_new/cell_topic.csv',
  'grn_extended': 'output/sp_new/grn_extended.csv'
}
## VIASH END
import argparse
argparser = argparse.ArgumentParser()
argparser.add_argument('--rna', type=str, help='Path to the input RNA data in h5ad format.')
argparser.add_argument('--atac', type=str, help='Path to the input ATAC data in h5ad format.')
argparser.add_argument('--prediction', type=str, help='Path to the output prediction in h5ad format.')
argparser.add_argument('--grn_extended', type=str, help='Path to the output prediction extended in csv format.')

args = argparser.parse_args()
if args.rna is not None:
    par['rna'] = args.rna
if args.atac is not None:
    par['atac'] = args.atac
if args.prediction is not None:
    par['prediction'] = args.prediction


try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
        'resources_dir': 'src/methods/multi_omics/scenicplus',
        'utils_dir': 'src/utils'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["utils_dir"])
from helper import download_databases, process_peak, run_cistopic, process_topics, preprocess_rna, snakemake_pipeline, post_process 

def print_memory_usage():
    import psutil

    process = psutil.Process(os.getpid())
    mem_info = process.memory_info().rss / (1024 * 1024)  # Convert to MB
    print(f"Memory usage: {mem_info:.2f} MB")

def main(par):
    # - set pathes
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
    par['chromsizes'] = f"{par['temp_dir']}/chromsizes.tsv"
    par['annotation_file'] = os.path.join(par['temp_dir'], 'gene_annotation.gtf')
    os.makedirs(par['atac_dir'], exist_ok=True)

    
    print('------- download_databases -------')
    download_databases(par)
    print_memory_usage()
    print('------- process_peak -------')
    process_peak(par)
    print_memory_usage()
    print('------- run_cistopic -------')
    run_cistopic(par)
    print_memory_usage()
    print('------- process_topics -------')
    process_topics(par)
    print_memory_usage()
    print('------- preprocess_rna -------')
    preprocess_rna(par)
    print_memory_usage()
    print('------- snakemake_pipeline -------')
    snakemake_pipeline(par)
    print_memory_usage()
    print('------- post_process -------')
    net = post_process(par)
    return net
if __name__ == '__main__':
    # - subset to one donor for test
    os.makedirs(par['temp_dir'], exist_ok=True)
    if False: #TODO: remove this
        adata = ad.read(par['rna'])
        adata = adata[adata.obs['donor_id']=='donor_0']
        adata.obs['donor_id'] = adata.obs['donor_id'].astype(str)
        assert adata.shape[0]>0, 'no cell left after filtering'
        par['rna'] = f"{par['temp_dir']}/rna.h5ad"
        adata.write(par['rna'])

        adata = ad.read(par['atac'])
        adata = adata[adata.obs['donor_id']=='donor_0']
        adata.obs['donor_id'] = adata.obs['donor_id'].astype(str)
        assert adata.shape[0]>0, 'no cell left after filtering'
        par['atac'] = f"{par['temp_dir']}/atac.h5ad"
        adata.write(par['atac'])
    # - main 
    net = main(par)
    
    dataset_id = ad.read_h5ad(par['rna'], backed='r').uns['dataset_id']
    output = ad.AnnData(X=None, uns={"method_id": 'scenicplus', "dataset_id": dataset_id, "prediction": net[["source", "target", "weight"]]})
    output.write(par['prediction'])
