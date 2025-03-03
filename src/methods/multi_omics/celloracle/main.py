import pandas as pd
import os 
from celloracle import motif_analysis as ma
import pandas as pd
import celloracle as co
import anndata
import numpy as np
import anndata as ad
from celloracle import motif_analysis as ma
import genomepy
import scanpy as sc

def base_grn(par) -> None:
    print("Reading atac data")
    atac = ad.read_h5ad(par["atac"])

    print("Format peak data")
    peaks = atac.var_names.to_numpy()
    peaks = [peak.replace(':','_').replace("-",'_') for peak in peaks]
    tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome="hg38")
    tss_annotated['peak_id'] = tss_annotated['chr'].astype(str)+"_"+tss_annotated['start'].astype(str)+"_"+tss_annotated['end'].astype(str)
    peak_gene = tss_annotated

    try:
        print("Install ref genome")
        genomepy.install_genome(name="hg38", provider="UCSC", genomes_dir=None)
    except:
        print("Couldnt install genome. Will look for the default location")

    ref_genome = "hg38"

    genome_installation = ma.is_genome_installed(ref_genome=ref_genome,
                                                genomes_dir=None)
    print(ref_genome, "installation: ", genome_installation)

    print("Instantiate TFinfo object")
    tfi = ma.TFinfo(peak_data_frame=peak_gene, 
                    ref_genome=ref_genome,
                    genomes_dir=None) 
    print("Motif scan")
    tfi.scan(fpr=0.05, 
            motifs=None,  # If you enter None, default motifs will be loaded.
            verbose=True)
    # Reset filtering 
    tfi.reset_filtering()

    print("Do filtering")
    tfi.filter_motifs_by_score(threshold=10)

    print("Format post-filtering results.")
    tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

    df = tfi.to_dataframe()
    print("Base GRN is built")
    df.to_csv(par['base_grn'])

def preprocess_rna(par) -> None:
    print("Processing rna data")
    adata = ad.read_h5ad(par['rna'])
    if True: #only one cluster
        adata.obs['cell_type'] = 'one_cell_type'
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')
    n_top_genes = min([3000, adata.shape[1]])

    filter_result = sc.pp.filter_genes_dispersion(adata.X,
                                                flavor='cell_ranger',
                                                n_top_genes=n_top_genes,
                                                log=False)

    # Subset the genes
    if False: #only hvgs
        adata = adata[:, filter_result.gene_subset]
    else:
        pass

    # Renormalize after filtering
    sc.pp.normalize_per_cell(adata)

    # Log transformation and scaling
    sc.pp.log1p(adata)
    sc.pp.scale(adata)

    # PCA
    sc.tl.pca(adata, svd_solver='arpack')

    # Diffusion map
    sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)

    sc.tl.diffmap(adata)
    # Calculate neihbors again based on diffusionmap 
    sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')

    sc.tl.louvain(adata, resolution=0.8)

    sc.tl.paga(adata, groups='louvain')

    sc.pl.paga(adata)
    sc.tl.draw_graph(adata, init_pos='paga', random_state=123)
    sc.pl.draw_graph(adata, color='louvain', legend_loc='on data')
    # Check data in anndata
    print("Metadata columns :", list(adata.obs.columns))
    print("Dimensional reduction: ", list(adata.obsm.keys()))

    adata.X = adata.layers["counts"].copy()

    adata.write_h5ad(f"{par['temp_dir']}/adata.h5ad")

def infer_grn(par):
    print("Inferring GRN using base GRN and rna expression")
    adata = ad.read_h5ad(f"{par['temp_dir']}/adata.h5ad")
    base_GRN = pd.read_csv(par['base_grn'])
    # Instantiate Oracle object
    oracle = co.Oracle()
    # Instantiate Oracle object.
    oracle.import_anndata_as_raw_count(adata=adata,
                                    cluster_column_name="cell_type",
                                    embedding_name="X_draw_graph_fr")
    # You can load TF info dataframe with the following code.
    oracle.import_TF_data(TF_info_matrix=base_GRN)

    oracle.perform_PCA()
    n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
    print(n_comps)
    n_comps = min(n_comps, 50)

    n_cell = oracle.adata.shape[0]
    print(f"cell number is :{n_cell}")
    k = min([int(0.025*n_cell), 50])
    print(f"Auto-selected k is :{k}")
    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                        b_maxl=k*4, n_jobs=par["num_workers"])
    links = oracle.get_links(cluster_name_for_GRN_unit="cell_type", alpha=10,
                        verbose_level=10, n_jobs=par["num_workers"])
    links.to_hdf5(file_path=par['links'])
def refine_grns(par):
    print("Refining GRNs")
    links_o = co.load_hdf5(par['links'])
    links_dict =  links_o.links_dict.copy()
    grn_stack = []
    tt = 0.05
    for cell_type, grn in links_dict.items():
        print(f"{cell_type}, GRN before filter: {grn.shape}")
        mask = grn.p<tt # remove insig links
        grn = grn[mask]
        grn = grn[~(grn.coef_abs==0)] # remove those with 0 coeff
        # filter based on z score 
        # z_scores = (grn.coef_abs - grn.coef_abs.mean())/grn.coef_abs.std()
        # mask = z_scores > 2
        # Sort by absolute coefficient values
        grn = grn.sort_values(by="coef_abs", ascending=False)

        # Select the top 50,000 links based on absolute weight
        mask = grn.index[:par['max_n_links']]
        grn = grn.loc[mask, :]

        grn = grn[['source', 'target', 'coef_mean']]

        grn.columns = ['source', 'target', 'weight']

        print(cell_type, 'links:', len(grn), ' TFs: ', grn.source.unique().shape[0], ' target: ', grn.target.unique().shape[0],)    
        grn['cell_type'] = cell_type
        grn_stack.append(grn)
    grn = pd.concat(grn_stack).reset_index(drop=True)
    grn = grn.groupby(['source', 'target'])['weight'].apply(np.mean).to_frame().reset_index()
    return grn

def main(par):
    os.makedirs(par['temp_dir'], exist_ok=True)
    base_grn(par)
    preprocess_rna(par)
    infer_grn(par)
    grn = refine_grns(par)
    return grn
    