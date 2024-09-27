
import pandas as pd
import anndata as ad
import sys
import json
import numpy as np
import matplotlib.pyplot as plt


def main(par):
    grn_model = pd.read_csv(par["grn_model"])
    
    exp_obj = Exp_analysis(grn_model)
    stats = exp_obj.calculate_basic_stats()
    
    if ('stats' in par) & (par['stats'] is not None):
        with open(par['stats'], 'w') as ff:
            json.dump(stats, ff)

    # print('Plotting CMD of weight: ', par['reg_weight_distribution'], flush=True)
    # fig, ax = plot_cumulative_density(grn_model.weight, title='CMD of reg. weight')
    # fig.savefig(par['reg_weight_distribution'])

    # print('Reading input files', flush=True)
    # # peak_gene_net['source'] = peak_gene_net['peak']
    # info_obj = Explanatory_analysis(net=grn_model)
    # print("Calculate basic stats")
    # stats = info_obj.calculate_basic_stats()
    # print(stats)
    # print("Outputting stats to :", par['stats'])
    
    # # print("Annotation of peaks")
    # # peak_annot = info_obj.annotate_peaks(annot_peak_database)
    # # print("Annotation of genes")
    # # gene_annot = info_obj.annotate_genes(annot_gene_database)
    # print("Topological analysis")
    # info_obj.calculate_centrality_stats()

    # tf_gene_in = info_obj.tf_gene.in_deg
    # tf_gene_out = info_obj.tf_gene.out_deg

    # print("Plotting tf-gene in degree, dir: ", par['tf_gene_indegee_fig'])
    # print("Plotting tf-gene out degree, dir: ", par['tf_gene_outdegee_fig'])
    # fig, ax = info_obj.plot_cdf(tf_gene_in, title='In degree TF-gene')
    # fig.savefig(par['tf_gene_indegee_fig'], dpi=300, bbox_inches='tight', format='png')
    # fig, ax = info_obj.plot_cdf(tf_gene_out, title='Out degree TF-gene')
    # fig.savefig(par['tf_gene_outdegee_fig'], dpi=300, bbox_inches='tight', format='png')