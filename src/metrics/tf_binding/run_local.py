import subprocess

import anndata as ad
import pandas as pd


for dataset in ["op"]:
    #for method in ['granie']:
    for method in ["negative_control", "granie", "ppcor", "portia", "pearson_corr", "positive_control", "scenic", "grnboost", "scprint", "scenicplus", "celloracle", "scglue", "figr"]:
        
        print()
        print(method)

        score_filepath = f"output/tf_binding/tf_binding.h5ad"
        subprocess.call([
            "python",
            "src/metrics/tf_binding/script.py",
            "--prediction", f"resources/results/{dataset}/{dataset}.{method}.{method}.prediction.h5ad",
            "--evaluation_data", f"resources/grn_benchmark/evaluation_data/{dataset}_bulk.h5ad",
            "--ground_truth", "resources/grn_benchmark/ground_truth/K562.csv",
            "--score", score_filepath
        ])

        adata = ad.read_h5ad(score_filepath)
        if "metric_values" in adata.uns:
            metric_names = adata.uns["metric_ids"]
            metric_values = adata.uns["metric_values"]
            df = pd.DataFrame({"metric": metric_names, "value": metric_values})
            df["dataset"] = dataset
            df["method"] = method
            df = df[["dataset", "method", "metric", "value"]]  # Reorder columns to match header
            df.to_csv(f"output/tf_binding/tf_binding_scores_{dataset}.csv", mode="a", header=False, index=False)
