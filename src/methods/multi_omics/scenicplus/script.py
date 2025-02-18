
import sys 
import os
import anndata as ad
## VIASH START
par = {
  'rna': 'resources_test/grn_benchmark/inference_data/op_rna.h5ad',
  'atac': 'resources_test/grn_benchmark/inference_data/op_atac.h5ad',
  'temp_dir': 'output/scenicplus',
  'prediction': 'output/prediction.h5ad',
  'qc': False,
  'num_workers': 20,
  'scplus_mdata': 'output/scenicplus_new/scplus_mdata.h5mu',
  'cell_topic': 'output/scenicplus_new/cell_topic.csv',
  'grn_extended': 'output/scenicplus_new/grn_extended.csv'
}
## VIASH END

import argparse
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--atac', type=str, help='Path to the multiomics atac file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=int, help='Number of cores')
parser.add_argument('--max_n_links', type=int)
args = parser.parse_args()

for key, value in vars(args).items():
    if value:
        par[key] = value

print(par)

try:
    sys.path.append(meta["resources_dir"])
except:
    pass
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
    net = main(par)
    dataset_id = ad.read_h5ad(par['rna'], backed='r').uns['dataset_id']
    output = ad.AnnData(X=None, uns={"method_id": 'scenicplus', "dataset_id": dataset_id, "prediction": net[["source", "target", "weight"]]})
    output.write(par['prediction'])
