import argparse
import sys

import anndata as ad
import numpy as np
import pandas as pd
import torch
from huggingface_hub import hf_hub_download
from scdataloader import Preprocessor
from scdataloader.utils import load_genes
from scipy.sparse import csr_matrix
from scprint import scPrint
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
    "how": "most var across",
}
## VIASH END
try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {"resources_dir": "src/utils"}
    sys.path.append(meta["resources_dir"])
from util import parse_args, process_links, efficient_melting
par = parse_args(par)

adata = ad.read_h5ad(par["rna"], backed="r")
train_perturbs = adata.obs['perturbation'].unique()
if adata.uns['dataset_id'] in ['replogle', 'xaira_HCT116', 'xaira_HEK293T']:
    tf_all = np.loadtxt(par['tf_all'], dtype=str)
    train_perturbs = np.intersect1d(tf_all, train_perturbs) 
    train_perturbs = train_perturbs[:100]  # limit to 100 perturbations
elif adata.uns['dataset_id'] in ['parsebioscience']:
    train_perturbs = train_perturbs[:10]
mask = adata.obs['perturbation'].isin(train_perturbs)
adata = adata[mask].to_memory()
    
adata.obs["is_primary_data"] = True
adata.obs["organism_ontology_term_id"] = "NCBITaxon:9606"
# adata.var = adata.var.set_index("gene_ids")

preprocessor = Preprocessor(
    min_valid_genes_id=min(0.7 * adata.n_vars, 10000),  # 90% of features up to 10,000
    # Turn off cell filtering to return results for all cells
    filter_cell_by_counts=400,
    min_nnz_genes=200,
    do_postp=False,
    is_symbol=True,
    force_preprocess=True,
    # Skip ontology checks
    skip_validate=True,
    use_raw=False,
)
if adata.raw is not None and adata.raw.X.shape[1] != adata.X.shape[1]:
    print("removing raw")
    del adata.raw
if adata.layers is not None and "counts" in adata.layers:
    adata.X = adata.layers["counts"]
    del adata.layers["counts"]


dataset_id = adata.uns["dataset_id"] if "dataset_id" in adata.uns else par["dataset_id"]

print("\n>>> Preprocessing data...", flush=True)
adata = preprocessor(adata)

if adata[0].X.sum() != int(adata[0].X.sum()):
    print("WARNING: you are not using count data")
    print("reverting logp1")
    adata.X = csr_matrix(np.power(adata.X.todense(), 2) - 1)


model_checkpoint_file = par["model"]
if model_checkpoint_file is None:
    print(f"\n>>> Downloading '{par['model_name']}' model...", flush=True)
    model_checkpoint_file = hf_hub_download(
        repo_id="jkobject/scPRINT", filename=f"{par['model_name']}.ckpt"
    )

if torch.cuda.is_available():
    print("CUDA is available, using GPU", flush=True)
    dtype = torch.float16
    transformer = "flash"
else:
    print("CUDA is not available, using CPU", flush=True)
    dtype = torch.float32
    transformer = "normal"

print(f"Model checkpoint file: '{model_checkpoint_file}'", flush=True)

# make sure that you check if you have a GPU with flashattention or not (see README)
try:
    m = torch.load(model_checkpoint_file)
# if not use this instead since the model weights are by default mapped to GPU types
except RuntimeError:
    m = torch.load(model_checkpoint_file, map_location=torch.device("cpu"))

if "prenorm" in m["hyper_parameters"]:
    m["hyper_parameters"].pop("prenorm")
    torch.save(m, model_checkpoint_file)

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

missing = set(model.genes) - set(load_genes(model.organisms).index)
if len(missing) > 0:
    print(
        "Warning: some genes missmatch exist between model and ontology: solving...",
    )
    model._rm_genes(missing)

if model.device.type == "cpu" and torch.cuda.is_available():
    print("CUDA is available, moving model to GPU", flush=True)
    model = model.to("cuda")


def main_sub(adata, model, par, cell_type):
    try:
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
            dtype=dtype,
        )

        grn = grn_inferer(model, adata, cell_type=cell_type)
    except ValueError as e:
        if "Extrapolation" in str(e):
            print(f"WARNING: {cell_type} has poor quality expression data")
            print("rerunning with most expressed genes only")
            grn_inferer.how = "most expr"
            grn = grn_inferer(model, adata, cell_type=cell_type)
        else:
            raise e
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


if "cell_type" not in adata.obs:
    adata.obs["cell_type"] = "dummy_cell_type"
    par["how"] = "most var within"
for i, cell_type in enumerate(adata.obs["cell_type"].unique()):
    print(cell_type)
    net = main_sub(adata, model, par, cell_type)
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
