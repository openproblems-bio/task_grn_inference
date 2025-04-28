import numpy as np
import pandas as pd
from arboreto.algo import grnboost2
from distributed import Client, LocalCluster
from tqdm import tqdm
import sys
import anndata as ad
import argparse

## VIASH START
par = {
    "rna": "resources/grn-benchmark/rna_d0_hvg.h5ad",
    "tf_all": "resources/grn_benchmark/prior/tf_all.csv",
    "prediction": "output/grnboost2_donor_0_hvg.csv",
    "max_n_links": 50000,
    "cell_type_specific": False,
    "normalize": False,
    "qc": False,
    "layer": "X_norm",
}
## VIASH END

## LOCAL START
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument("--rna", type=str, help="Path to the multiomics RNA file")
parser.add_argument("--prediction", type=str, help="Path to the prediction file")
parser.add_argument("--resources_dir", type=str, help="Path to the prediction file")
parser.add_argument("--tf_all", type=str, help="Path to the tf_all")
parser.add_argument("--max_n_links", type=str, help="Number of top links to retain")
parser.add_argument("--dataset_id", type=str, help="Dataset id")
parser.add_argument("--normalize", action="store_true")
args = parser.parse_args()

par_local = vars(args)

for key, value in par_local.items():
    if value is not None:
        par[key] = value

## LOCAL END

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {"resources_dir": "src/utils"}
    sys.path.append(meta["resources_dir"])

from util import process_links, basic_qc


def main(par: dict) -> pd.DataFrame:
    """
    Main function to infer GRN
    """
    print("Reading data")
    adata_rna = ad.read_h5ad(par["rna"])
    if "qc" in par:
        if par["qc"]:
            print("Shape before QC: ", adata_rna.shape)
            adata_rna = basic_qc(adata_rna)
            print("Shape after QC: ", adata_rna.shape)

    gene_names = adata_rna.var_names

    X = adata_rna.layers[par["layer"]]

    # Load list of putative TFs
    tfs = np.loadtxt(par["tf_all"], dtype=str)
    tf_names = [gene_name for gene_name in gene_names if (gene_name in tfs)]

    # GRN inference
    client = Client(processes=False)

    print("Infer grn", flush=True)

    network = grnboost2(
        X, client_or_address=client, gene_names=gene_names, tf_names=tf_names
    )
    network.rename(
        columns={"TF": "source", "target": "target", "importance": "weight"},
        inplace=True,
    )
    network.reset_index(drop=True, inplace=True)
    network = process_links(network, par)
    return network


if __name__ == "__main__":
    dataset_id = ad.read_h5ad(par["rna"], backed="r").uns["dataset_id"]

    net = main(par)

    # Save inferred GRN
    print("Output GRN")
    # convert the predictions to the benchmark format
    net = net.astype(str)
    output = ad.AnnData(
        X=None,
        uns={
            "method_id": "grnboost2",
            "dataset_id": dataset_id,
            "prediction": net[["source", "target", "weight"]],
        },
    )
    output.write(par["prediction"])
