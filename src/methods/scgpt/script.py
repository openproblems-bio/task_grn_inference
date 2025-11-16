import os
import sys

import anndata as ad
import gdown
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scipy.sparse import csr_matrix

## VIASH START
par = {
    "rna": "resources_test/grn_benchmark/inference_data/op_rna.h5ad",
    "tf_all": "resources_test/prior/tf_all.csv",
    "prediction": "output/scgpt/prediction.h5ad",
    "max_n_links": 50000,
    "batch_size": 16,
    "temp_dir": "output/scgpt/",
    "file_dir": "src/methods/scgpt",
    "num_genes": 3000,
    "max_cells": 1000,
    "num_attn_layers": 11,
}
## VIASH END
try:
    sys.path.append(meta["resources_dir"])
    print(meta["resources_dir"])
except:
    meta = {"resources_dir": "src/utils"}
    sys.path.append(meta["resources_dir"])
from scgpt_helper import generate_grn, prepare_model
from util import efficient_melting, parse_args

par = parse_args(par)

## GETTING DATA

os.environ["KMP_WARNINGS"] = "off"


if __name__ == "__main__":

    adata = sc.read(par["rna"], backed="r")

    if adata.uns["dataset_id"] in ["replogle", "xaira_HCT116", "xaira_HEK293T"]:
        train_perturbs = adata.obs["perturbation"].unique()
        tf_all = np.loadtxt(par["tf_all"], dtype=str)
        train_perturbs = np.intersect1d(tf_all, train_perturbs)
        train_perturbs = train_perturbs[:100]  # limit to 100 perturbations
        mask = adata.obs["perturbation"].isin(train_perturbs)
        adata = adata[mask].to_memory()
    elif adata.uns["dataset_id"] in ["parsebioscience"]:
        train_perturbs = adata.obs["perturbation"].unique()
        train_perturbs = train_perturbs[:10]
        mask = adata.obs["perturbation"].isin(train_perturbs)
        adata = adata[mask].to_memory()
    else:
        adata = adata.to_memory()

    if adata.raw is not None and adata.raw.X.shape[1] != adata.X.shape[1]:
        print("removing raw")
        del adata.raw
    if adata.layers is not None and "counts" in adata.layers:
        adata.X = adata.layers["counts"]
        del adata.layers["counts"]

    if adata[0].X.sum() != int(adata[0].X.sum()):
        print("WARNING: you are not using count data")
        print("reverting logp1")
        from scipy.sparse import issparse
        if issparse(adata.X):
            adata.X = csr_matrix(np.expm1(adata.X.toarray()))
        else:
            adata.X = csr_matrix(np.expm1(adata.X))

    adata.var["symbol"] = adata.var.index
    if "gene_id" not in adata.var.columns:
        from util import add_gene_id
        adata = add_gene_id(adata)
    adata.var["ensembl_id"] = adata.var["gene_id"].values
    dataset_id = adata.uns["dataset_id"] if "dataset_id" in adata.uns else par["dataset_id"]


    ### LOADING MODEL AND PREDICTION
    os.makedirs(par["temp_dir"], exist_ok=True)
    # Download datasets
    model_file = f"{par['temp_dir']}/best_model.pt"

    print(os.system("pip list | grep flas"))

    if not os.path.exists(model_file):
        gdown.download(
            "https://drive.google.com/uc?id=1MJaavaG0ZZkC_yPO4giGRnuCe3F1zt30",
            par["temp_dir"] if par["temp_dir"].endswith("/") else par["temp_dir"] + "/",
            quiet=False,
        )
        gdown.download(
            "https://drive.google.com/uc?id=127FdcUyY1EM7rQfAS0YI4ms6LwjmnT9J",
            par["temp_dir"] if par["temp_dir"].endswith("/") else par["temp_dir"] + "/",
            quiet=False,
        )
        gdown.download(
            "https://drive.google.com/uc?id=1y4UJVflGl-b2qm-fvpxIoQ3XcC2umjj0",
            par["temp_dir"] if par["temp_dir"].endswith("/") else par["temp_dir"] + "/",
            quiet=False,
        )

    print(os.listdir(par["temp_dir"]))
    model, vocab = prepare_model(model_dir=par["temp_dir"])

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)

    #### PREDICTION

    if "cell_type" not in adata.obs:
        adata.obs["cell_type"] = "dummy_cell_type"
        par["how"] = "most var within"
    for i, cell_type in enumerate(adata.obs["cell_type"].unique()):
        print(cell_type)
        subadata = adata[adata.obs["cell_type"] == cell_type]
        # Handle both sparse and dense matrices
        gene_sums = subadata.X.sum(0)
        if hasattr(gene_sums, 'A1'):
            gene_sums = gene_sums.A1
        else:
            gene_sums = np.asarray(gene_sums).ravel()
        subadata = subadata[
            : par["max_cells"], gene_sums.argsort()[::-1][: par["num_genes"]]
        ]
        subadata, net = generate_grn(
            model,
            vocab,
            subadata,
            batch_size=par["batch_size"],
            num_attn_layers=par["num_attn_layers"],
        )
        gene_names = subadata.var["symbol"].values
        net = efficient_melting(net.T, gene_names, symmetric=False)
        net = net[net["weight"] != 0]

        # - subset to TFs
        tfs = set(np.loadtxt(par["tf_all"], dtype=str))

        net = net[net["source"].isin(tfs)]

        net = net.sort_values(
            by="weight", ascending=False, key=abs
        )[
            : 2 * par["max_n_links"]
        ]  # i set this to double of allowed link just in case the symmetry exists. metrics will take care of this
        net["cell_type"] = cell_type
        if i == 0:
            net_all = net
        else:
            net_all = pd.concat([net_all, net], axis=0)


    print(f"Writing results to {par['prediction']}")
    net_all = net_all.astype(str)
    output = ad.AnnData(
        X=None,
        uns={
            "method_id": "scgpt",
            "dataset_id": dataset_id,
            "prediction": net_all[["source", "target", "weight", "cell_type"]],
        },
    )
    output.write(par["prediction"])
