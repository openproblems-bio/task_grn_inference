import sys
from huggingface_hub import hf_hub_download
import numpy as np
import anndata as ad
import pandas as pd
import torch
import argparse
from scprint import scPrint
from scdataloader import Preprocessor
from scprint.tasks import GNInfer

torch.set_float32_matmul_precision("medium")

## VIASH START
par = {
    "rna": "resources/grn_benchmark/inference_data/op_rna.h5ad",
    "tf_all": "resources/grn_benchmark/prior/tf_all.csv",
    "prediction": "output/grn.h5ad",
    "filtration": "none",  #'top-k',
    "max_n_links": 50_000,
    "num_genes": 5000,
    "max_cells": 1000,
    "num_workers": 20,
    "model_name": "v2-medium",
    "model": None,
    "dataset_id": "op",
    "is_test": True,
    "how": "most var within",
}
## VIASH END

## LOCAL START
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument("--rna", type=str, help="Path to the multiomics RNA file")
parser.add_argument("--prediction", type=str, help="Path to the prediction file")
parser.add_argument("--resources_dir", type=str, help="Path to the prediction file")
parser.add_argument("--tf_all", type=str, help="Path to the tf_all")
parser.add_argument("--num_workers", type=int, help="Number of cores")
parser.add_argument("--max_n_links", type=int, help="Number of top links to retain")
parser.add_argument("--dataset_id", type=str, help="Dataset id")
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
from util import efficient_melting


def main_sub(adata, model, par):
    grn_inferer = GNInfer(
        how=par["how"],
        preprocess="softmax",
        head_agg="mean",
        forward_mode="none",
        filtration=par["filtration"],
        num_genes=par["num_genes"],
        num_workers=par["num_workers"],
        max_cells=par["max_cells"],
        doplot=False,
        batch_size=16,
    )

    grn = grn_inferer(model, adata)

    net = grn.varp["GRN"]
    if hasattr(net, "todense"):  # Check if it's a sparse matrix
        net = net.todense().A.T
    else:  # Already a NumPy array
        net = net.T

    # - melt to have the format of the benchmark
    gene_names = grn.var["symbol"].values
    net = efficient_melting(net, gene_names, symmetric=False)
    net = net[net["weight"] != 0]

    # - subset to TFs
    tfs = set(np.loadtxt(par["tf_all"], dtype=str))

    net = net[net["source"].isin(tfs)]

    net = net.sort_values(
        by="weight", ascending=False, key=abs
    )[
        : 2 * par["max_n_links"]
    ]  # i set this to double of allowed link just in case the symmetry exists. metrics will take care of this

    return net


adata = ad.read_h5ad(par["rna"])
adata.obs["is_primary_data"] = True
adata.obs["organism_ontology_term_id"] = "NCBITaxon:9606"
# adata.var = adata.var.set_index("gene_ids")

print("\n>>> Preprocessing data...", flush=True)
preprocessor = Preprocessor(
    min_valid_genes_id=min(0.7 * adata.n_vars, 10000),  # 90% of features up to 10,000
    # Turn off cell filtering to return results for all cells
    filter_cell_by_counts=400,
    min_nnz_genes=200,
    do_postp=False,
    is_symbol=True,
    # Skip ontology checks
    skip_validate=True,
)

dataset_id = adata.uns["dataset_id"] if "dataset_id" in adata.uns else par["dataset_id"]
adata = preprocessor(adata)

model_checkpoint_file = par["model"]
if model_checkpoint_file is None:
    print(f"\n>>> Downloading '{par['model_name']}' model...", flush=True)
    model_checkpoint_file = hf_hub_download(
        repo_id="jkobject/scPRINT", filename=f"{par['model_name']}.ckpt"
    )

if torch.cuda.is_available():
    print("CUDA is available, using GPU", flush=True)
    precision = "16"
    dtype = torch.float16
    transformer = "flash"
else:
    print("CUDA is not available, using CPU", flush=True)
    precision = "32"
    dtype = torch.float32
    transformer = "normal"

print(f"Model checkpoint file: '{model_checkpoint_file}'", flush=True)

m = torch.load(model_checkpoint_file, map_location=torch.device("cpu"))
if "label_counts" in m["hyper_parameters"]:
    model = scPrint.load_from_checkpoint(
        model_checkpoint_file,
        transformer=transformer,  # Don't use this for GPUs with flashattention
        precpt_gene_emb=None,
        classes=m["hyper_parameters"]["label_counts"],
    )
else:
    model = scPrint.load_from_checkpoint(
        model_checkpoint_file,
        transformer=transformer,  # Don't use this for GPUs with flashattention
        precpt_gene_emb=None,
    )
del m

if model.device.type == "cpu" and torch.cuda.is_available():
    print("CUDA is available, moving model to GPU", flush=True)
    model = model.to("cuda")

if "cell_type" not in adata.obs:
    adata.obs["cell_type"] = "dummy_cell_type"
for i, cell_type in enumerate(adata.obs["cell_type"].unique()):
    print(cell_type)
    adata_cell_type = adata[adata.obs["cell_type"] == cell_type].copy()
    net = main_sub(adata_cell_type, model, par)
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
        "method_id": "scprint",
        "dataset_id": dataset_id,
        "prediction": net_all[["source", "target", "weight", "cell_type"]],
    },
)
output.write(par["prediction"])
