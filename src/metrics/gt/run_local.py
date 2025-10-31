import subprocess

import anndata as ad
import pandas as pd


results = {"dataset": [], "method": [], "score": []}
for dataset in ["beeline_hESC", "beeline_hHep", "beeline_mDC", "beeline_mESC", "beeline_mESC-E", "beeline_mHSC-E", "beeline_mHSC-L"]:
    #for method in ["negative_control", "granie", "ppcor", "portia", "pearson_corr", "positive_control", "scenic", "grnboost", "scprint", "scenicplus", "celloracle", "scglue", "figr"]:
    #for method in ['granie', 'ppcor', 'portia', 'negative_control']:
    for method in ['portia']:
        
        print()
        print(method)

        score_filepath = f"output/gt/gt.h5ad"
        if "_h" in dataset:
            gt = f"human/{dataset.split('_')[1]}-ChIP-seq-network"
        else:
            gt = f"mouse/{dataset.split('_')[1]}-ChIP-seq-network"
        subprocess.call([
            "python",
            "src/metrics/gt/script.py",
            "--prediction", f"resources/results/{dataset}/{method}.h5ad",
            "--evaluation_data", f'resources/grn_benchmark/inference_data/{dataset}.h5ad',
            "--ground_truth", f"resources/grn_benchmark/ground_truth/beeline/Networks/{gt}.csv",
            "--score", score_filepath
        ])

        adata = ad.read_h5ad(score_filepath)
        if "metric_values" in adata.uns:
            metric_names = adata.uns["metric_ids"]
            metric_values = adata.uns["metric_values"]
            print(metric_values)

df = pd.DataFrame(results)
df.to_csv(f"output/gt/gt_scores_beeline.csv", header=True, index=False)
