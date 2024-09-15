import sctk
import pandas as pd 
import anndata as ad 
import scanpy as sc
import matplotlib.pyplot as plt
import os
par = {
    "multiomics_rna":"resources/grn-benchmark/multiomics_rna.h5ad"
}
multiomics_rna = ad.read_h5ad(par['multiomics_rna'])
# sctk.calculate_qc(multiomics_rna, flags={"mito": r"^MT-", "ribo": r"^RP[LS]"})
sctk.calculate_qc(multiomics_rna)
sctk.cellwise_qc(multiomics_rna)

os.makedirs('output/qc', exist_ok=True)
sc.set_figure_params(figsize=(6,6))

# plot n counts vs n genes
save_name = "output/qc/scatter_plot.png"
print("Plotting n-counts vs n-genes to")
p3 = sc.pl.scatter(multiomics_rna, "n_counts", "n_genes", 
                   color="cell_type", size=20, alpha=.5)

# min max of different metrics
def min_max(name, data):
    print(f'{name}: min:{data.min()}, max: {data.max()}')
min_max('n_counts', multiomics_rna.obs.n_counts)
min_max('n_genes', multiomics_rna.obs.n_genes)
min_max('percent_hb', multiomics_rna.obs.percent_hb)
min_max('percent_mito', multiomics_rna.obs.percent_mito)

# Plot CDF
save_name = "output/qc/cdf.png"
print("Saving cdf graphs to :", save_name)
def plot_CDF(data, title, ax):
    x = np.sort(data)
    y = np.arange(1, len(data) + 1) / len(data)
    ax.plot(x, y, marker='.', linestyle='none')
    ax.set_xlabel('Value')
    ax.set_ylabel('Cumulative Probability')
    ax.set_title(title)
    ax.grid(True)
fig, axes = plt.subplots(1, 5, figsize=(20,4))
plot_CDF(multiomics_rna.obs.n_counts, 'n_counts', axes[0])
plot_CDF(multiomics_rna.obs.n_genes, 'n_genes', axes[1])
plot_CDF(multiomics_rna.obs.percent_hb, 'percent_hb', axes[2])
plot_CDF(multiomics_rna.obs.percent_mito, 'percent_mito', axes[3])
plot_CDF(multiomics_rna.obs.percent_ribo, 'percent_ribo', axes[4])
plt.tight_layout()
fig.savefig(save_name, dpi=300)
