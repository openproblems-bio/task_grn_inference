import pandas as pd
import anndata as ad
import sys

## VIASH START
par = {
  "perturbation_data": "resources/grn-benchmark/pertubation_data.h5ad",
  "tf_gene_net": "resources/grn-benchmark/grn_models/figr.csv",
  "peak_gene_net": "resources/grn-benchmark/peak_gene_models/figr.csv",
  "annot_peak_database": "resources/grn-benchmark/supp/annot_peak_database.csv",
  "annot_gene_database": "resources/grn-benchmark/supp/annot_gene_database.csv",
  "hvgs":  "resources/grn-benchmark/supp/hvgs.txt"

}
## VIASH END

# meta = {
#   "resources_dir":'resources'
# }
sys.path.append(meta["resources_dir"])
from helper import Explanatory_analysis

print('Reading input files', flush=True)

perturbation_data = pd.read_csv(par["perturbation_data"])
tf_gene_net = pd.read_csv(par["tf_gene_net"])
peak_gene_net = pd.read_csv(par["peak_gene_net"])
annot_peak_database = pd.read_csv(par["annot_peak_database"])
hvgs = pd.read_csv(par["hvgs"])

print("Create exp anal object")
peak_gene_net['source'] = peak_gene_net['peak']
info_obj = Explanatory_analysis(net=tf_gene_net, peak_gene_net=peak_gene_net)
print("Annotation of peaks")
peak_annot = info_obj.annotate_peaks(annot_peak_database)
# print("Annotation of genes")
gene_annot = info_obj.annotate_genes(annot_gene_database)
print("Topological analysis")
info_obj.calculate_centrality_stats()
peak_gene = info_obj.peak_gene.in_deg
tf_gene_in = info_obj.tf_gene.in_deg
tf_gene_out = info_obj.tf_gene.out_deg

print("HVGs analysis")
pert_genes = perturbation_data.var_names
hvgs_info = info_obj.calculate_hvgs_stats(hvgs, pert_genes)



