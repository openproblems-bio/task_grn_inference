import pandas as pd
import anndata as ad
import sys
import os
import requests
## VIASH START
par = {
  "multiomics_rna": "resources/grn-benchmark/multiomics_rna_0.h5ad",
  "multiomics_atac": "resources/grn-benchmark/multiomics_atac_0.h5ad",
  "temp_dir": 'output/scglue/',
  "num_workers": 20,
  "prediction": "output/scglue/prediction.csv",
  "nes_threshold": 1.5,
  "rank_threshold": 1500,
  "top_n_targets": 100
}
## VIASH END
# meta = {
#   "resources_dir": 'src/methods/multi_omics/scglue/'
# }

def download_annotation(par):
    # get gene annotation
    response = requests.get("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz")
    par['annotation_file'] = f"{par['temp_dir']}/gencode.v45.annotation.gtf.gz"
    if response.status_code == 200:
        with open(par['annotation_file'], 'wb') as file:
            file.write(response.content)
        print(f"File downloaded and saved to {par['annotation_file']}")
    else:
        print(f"Failed to download the gencode.v45.annotation.gtf.gz. Status code: {response.status_code}")
def download_motifs(par):
    # get gene annotation
    tag = "JASPAR2022-hg38.bed.gz"
    response = requests.get(f"http://download.gao-lab.org/GLUE/cisreg/{tag}")
    par['motif_file'] = f"{par['temp_dir']}/{tag}"
    if response.status_code == 200:
        with open(par['motif_file'], 'wb') as file:
            file.write(response.content)
        print(f"File downloaded and saved to {par['motif_file']}")
    else:
        print(f"Failed to download the {tag}. Status code: {response.status_code}")
      
print("Downloading prior started")
download_annotation(par)
print("Downloading prior ended")

print("Downloading motifs started")
download_motifs(par)
print("Downloading motifs ended")

sys.path.append(meta["resources_dir"])
from main import main 
os.makedirs(par['temp_dir'], exist_ok=True)
prediction = main(par)

print('Write output to file', flush=True)
prediction.to_csv(par["prediction"])



