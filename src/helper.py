import yaml
import os
import pandas as pd
import anndata as ad 
import numpy as np
import scanpy as sc 
import sys 
import subprocess
import io
import warnings

def process_seqera_trace():
    model = 'celloracle'
    base_dir = f'resources/results/{model}'
    trace = pd.read_csv(f'{base_dir}/trace.txt', sep='\t')
    trace_seqera = process_trace_seqera(trace)
    trace_seqera.index = [model]

    map_dict = {'peak_rss':'MaxRSS', 'duration':'Elapsed'}

    trace_seqera = trace_seqera[map_dict.keys()]
    trace_seqera.columns = trace_seqera.columns.map(map_dict)

    trace_seqera

    # merge local and cluster dfs 
    map_dict = {'peak_rss':'MaxRSS', 'duration':'Elapsed'}

    df_local = df_local[map_dict.values()]

    df_res = pd.concat([trace_seqera, df_local], axis=0)

    df_res.columns = ['Peak memory (GB)', 'Duration (hour)']

    df_res
    

def analyse_meta_cells(task_grn_inference_dir):
    dataset = 'op' #'nakatake' #op', norman


    par = {
        'rna': f'{task_grn_inference_dir}/resources/grn_benchmark/inference_datasets/{dataset}_rna.h5ad',
        "evaluation_data": f"{task_grn_inference_dir}/resources/grn_benchmark/evaluation_datasets//{dataset}_perturbation.h5ad",

        'layer': 'X_norm',
        'consensus':  f'{task_grn_inference_dir}/resources/grn_benchmark/prior/{dataset}_consensus-num-regulators.json',
        'tf_all': f'{task_grn_inference_dir}/resources/grn_benchmark/prior/tf_all.csv',
        'apply_tf': True,
        'apply_skeleton': False,
        'verbose': 2,
        'max_n_links': 50_000,
        'temp_dir': 'output/metacells/',
        'subsample': -1,
        'num_workers': 4,
        'binarize': True,
        'reg_type': 'ridge'
    }

    os.makedirs(par['temp_dir'], exist_ok=True)
    # - imports 
    sys.path.append(f'{task_grn_inference_dir}/')
    sys.path.append(f'{task_grn_inference_dir}/src/utils')

    from src.metrics.regression_1.main import main as main_reg1
    from src.metrics.regression_2.main import main as main_reg2

    from util import corr_net

    # - read inputs and cluster with differen resolutions
    rna = ad.read_h5ad(par['rna'])
    rna.X = rna.layers[par['layer']]

    sc.pp.pca(rna)
    sc.pp.neighbors(rna)

    for res in range(1, 20, 2):
        sc.tl.leiden(rna, resolution=res, key_added=f'leiden_{res}')
    rna.write(f"{par['temp_dir']}/rna.h5ad")
    # - pseudobulkd and run per res
    gene_names = rna.var_names
    for i_run, res in enumerate([-1]+list(range(1, 20, 2))): 
        if res == -1:
            X = rna.X
        else:
            cluster_key = f'leiden_{res}'
            expression_data = rna.to_df()  # Converts AnnData object to a DataFrame (cells × genes)

            # Add the cluster assignments to the expression DataFrame
            expression_data[cluster_key] = rna.obs[cluster_key].values

            # Group by cluster and calculate the mean expression per gene
            pseudobulk = expression_data.groupby(cluster_key).mean()
            #- create anndata 
            adata_pseudobulk = ad.AnnData(
                X=pseudobulk.values,  # Gene expression matrix (clusters × genes)
                obs=pd.DataFrame(index=pseudobulk.index),  # Observations (clusters)
                var=pd.DataFrame(index=pseudobulk.columns)  # Variables (genes)
            )
            X = adata_pseudobulk.X
        # - run corr 
        net = corr_net(X, gene_names, par)
        par['prediction'] = f"{par['temp_dir']}/net_{res}.csv"
        net.to_csv(par['prediction'])

        # - regression 1 and 2
        scores_reg1 = main_reg1(par)

        scores_reg2 = main_reg2(par)

        scores = pd.concat([scores_reg1, scores_reg2], axis=1)
        scores.index = [res]

        if i_run == 0:
            scores_all = scores
        else:
            scores_all = pd.concat([scores_all, scores], axis=0)
        print(scores_all)
        scores_all.to_csv(f"{par['temp_dir']}/scores_all.csv")

def analyse_imputation(task_grn_inference_dir):
    dataset = 'op' #'nakatake' #op', norman


    par = {
        'rna': f'{task_grn_inference_dir}/resources/grn_benchmark/inference_datasets/{dataset}_rna.h5ad',
        "evaluation_data": f"{task_grn_inference_dir}/resources/grn_benchmark/evaluation_datasets//{dataset}_perturbation.h5ad",

        'layer': 'X_norm',
        'consensus':  f'{task_grn_inference_dir}/resources/grn_benchmark/prior/{dataset}_consensus-num-regulators.json',
        'tf_all': f'{task_grn_inference_dir}/resources/grn_benchmark/prior/tf_all.csv',
        'apply_tf': True,
        'apply_skeleton': False,
        'verbose': 2,
        'max_n_links': 50_000,
        'temp_dir': 'output/imputation/',
        'subsample': -1,
        'num_workers': 4,
        'binarize': True,
        'reg_type': 'ridge'
    }

    os.makedirs(par['temp_dir'], exist_ok=True)
    # - imports 
    sys.path.append(f'{task_grn_inference_dir}/')
    sys.path.append(f'{task_grn_inference_dir}/src/utils')

    from src.metrics.regression_1.main import main as main_reg1
    from src.metrics.regression_2.main import main as main_reg2

    from util import corr_net

    # - read inputs and cluster with differen resolutions
    rna = ad.read_h5ad(par['rna'])
    rna.X = rna.layers[par['layer']]

    sc.pp.pca(rna)
    sc.pp.neighbors(rna)

    X_original = rna.X.todense().A
    

    # - pseudobulkd and run per res
    gene_names = rna.var_names
    for i_run, method in enumerate(['magic', 'knn',  -1]): 
        if method == -1:
            X = rna.X.todense().A
        elif method == 'knn':
            X_original[X_original == 0] = np.nan
            from sklearn.impute import KNNImputer
            knn_imputer = KNNImputer(n_neighbors=10) 
            X = knn_imputer.fit_transform(X_original)
        elif method == 'softimpute':
            from fancyimpute import SoftImpute
            X = SoftImpute().fit_transform(X_original)
        elif method == 'magic':
            import scprep
            import magic

            magic_operator = magic.MAGIC()
            X = magic_operator.fit_transform(X_original)
        else:
            raise ValueError('define first')
        # - run corr 
        net = corr_net(X, gene_names, par)
        par['prediction'] = f"{par['temp_dir']}/net_{method}.csv"
        net.to_csv(par['prediction'])

        # - regression 1 and 2
        scores_reg1 = main_reg1(par)

        scores_reg2 = main_reg2(par)

        scores = pd.concat([scores_reg1, scores_reg2], axis=1)
        scores.index = [method]

        if i_run == 0:
            scores_all = scores
        else:
            scores_all = pd.concat([scores_all, scores], axis=0)
        print(scores_all)
        scores_all.to_csv(f"{par['temp_dir']}/scores_all.csv")

def analyse_corr_vs_tfmasked_corr(task_grn_inference_dir):
    for i_run, dataset in enumerate(['op', 'replogle2', 'nakatake', 'norman', 'adamson']):
        par = {
            'rna': f'{task_grn_inference_dir}/resources/grn_benchmark/inference_datasets/{dataset}_rna.h5ad',
            "evaluation_data": f"{task_grn_inference_dir}/resources/grn_benchmark/evaluation_datasets//{dataset}_perturbation.h5ad",

            'layer': 'X_norm',
            'consensus':  f'{task_grn_inference_dir}/resources/grn_benchmark/prior/{dataset}_consensus-num-regulators.json',
            'tf_all': f'{task_grn_inference_dir}/resources/grn_benchmark/prior/tf_all.csv',
            'apply_skeleton': False,
            'verbose': 2,
            'max_n_links': 50_000,
            'temp_dir': 'output/causal_vs_corr/',
            'subsample': -1,
            'num_workers': 4,
            'binarize': True,
            'reg_type': 'ridge'
        }

        os.makedirs(par['temp_dir'], exist_ok=True)
        # - imports 
        sys.path.append(f'{task_grn_inference_dir}/')
        sys.path.append(f'{task_grn_inference_dir}/src/utils')

        from src.metrics.regression_1.main import main as main_reg1
        from src.metrics.regression_2.main import main as main_reg2

        from util import efficient_melting, process_links

        # - read inputs and cluster with differen resolutions
        rna = ad.read_h5ad(par['rna'])
        gene_names = rna.var_names.to_numpy()
        try:
            rna.X = rna.layers[par['layer']]
        except:
            raise ValueError(f'{dataset} doesnt have X_norm')
        X = rna.X
        import scipy.sparse
        # Assert that rna.X is not a sparse matrix
        if scipy.sparse.issparse(X):
            X = X.todense().A

        def calculate_corr(X, par):
            net = np.corrcoef(X.T)
            
            # # Convert to a DataFrame with gene names as both row and column indices
            net = pd.DataFrame(net, index=gene_names, columns=gene_names)
            
            net = net.values

            net = efficient_melting(net, gene_names, par)
            if par['apply_tf']:
                tf_all = np.loadtxt(par['tf_all'], dtype=str)
                tf_all = np.intersect1d(tf_all, gene_names)
                print('TF subsetting')
                net = net[net.source.isin(tf_all)]
            
            net = process_links(net, par)
            return net
        def calculate_scores(par):
            scores_reg1 = main_reg1(par)
            scores_reg2 = main_reg2(par)
            scores = pd.concat([scores_reg1, scores_reg2], axis=1)
            return scores
        # - causal net 
        par['apply_tf'] = True
        net = calculate_corr(X, par)
        par['prediction'] = f"{par['temp_dir']}/net_causal.csv"
        net.to_csv(par['prediction'])
        scores_causal = calculate_scores(par)
        scores_causal['type'] = ['causal']
        # - corr net 
        par['apply_tf'] = False
        net = calculate_corr(X, par)
        par['prediction'] = f"{par['temp_dir']}/net_corr.csv"
        net.to_csv(par['prediction'])
        scores_corr = calculate_scores(par)
        scores_corr['type']= ['corr']

        scores_all = pd.concat([scores_causal, scores_corr], axis=0)
        
        scores_all['dataset'] = dataset
        if i_run == 0:
            scores_all_datasets = scores_all
        else:
            scores_all_datasets = pd.concat([scores_all, scores_all_datasets], axis=0)
        print(scores_all)
        scores_all_datasets.to_csv(f"{par['temp_dir']}/scores_all.csv")
def check_scores(par):
    df_scores = pd.read_csv(f"resources/scores/scgen_pearson-ridge.csv", index_col=0)
    df_all_n = (df_scores-df_scores.min(axis=0))/(df_scores.max(axis=0)-df_scores.min(axis=0))
    df_scores['rank'] = df_all_n.mean(axis=1).rank(ascending=False).astype(int)
    return df_scores

def calculate_scores():
    command = [
        "sbatch", 
        "scripts/sbatch/calculate_scores.sh"
    ] 
    # Print command to verify
    subprocess.run(command)
    

def run_consensus(par):
    model_names = ' '.join(par['methods'])
    command = [
        "python", 
        "src/metrics/regression_2/consensus/script.py", 
        "--models_dir", par['models_dir'], 
        "--models"
    ] + model_names.split()

    # Print command to verify
    subprocess.run(command)

def run_grn_seqera():
    # if False: # submit the job
    # !bash src/methods/multi_omics/celloracle/run.sh
    # if False: # get the results
    #     !aws s3 sync s3://openproblems-data/resources/grn/results/celloracle resources/results/celloracle
    # if False: # process celloracle
    #     # - prediction 
    #     # !mv resources/results/celloracle/celloracle.grn_inference_celloracle.prediction.csv resources/results/celloracle/celloracle.csv 
    #     !cp resources/results/celloracle/celloracle.csv  {par['models_dir']}/celloracle.csv
    #     !ls {par['models_dir']}
        
    #     # -  peak - gene
    #     # !mkdir resources/grn_models/peak_gene/
    #     df = pd.read_csv("resources/results/celloracle/output/celloracle/base_grn.csv")[['peak_id','gene_short_name']]
    #     df.columns = ['peak','target']
    #     df.to_csv('resources/grn_models/peak_gene/celloracle.csv')
    raise ValueError('define first')


def create_skeleton():
    skeleton_encode = pd.read_csv('output/skeleton_encode/tf2gene.csv', index_col=0)
    skeleton_jaspar = pd.read_csv('output/skeleton_jaspar/tf2gene.csv', index_col=0)
    skeleton = pd.concat([skeleton_encode, skeleton_jaspar], axis=0)
    skeleton = skeleton.drop_duplicates().reset_index(drop=True)
    skeleton.nunique()
    all_links = skeleton['source'].astype(str) + '_' + skeleton['target'].astype(str)
    np.savetxt('resources/grn_benchmark/prior/skeleton.csv', all_links.values, fmt='%s')

def merge_tf_motifs():
    cols = ['chrom', 'chromStart',  'chromEnd', 'name']
    if True:
        
        # merge two db
        df_jaspar = pd.read_csv('output/db/JASPAR2022-hg38.bed.gz', sep='\t', names=cols, comment='#')

        df_encode = pd.read_csv('output/db/ENCODE-TF-ChIP-hg38.bed.gz', sep='\t', names=cols, comment='#')
        if False: #exp analysis
            df_jaspar['peaks'] = df_jaspar['chr']+df_jaspar['start'].astype(str)+df_jaspar['end'].astype(str)
            df_encode['peaks'] = df_encode['chr']+df_encode['start'].astype(str)+df_encode['end'].astype(str)

            print(df_jaspar['tf'].nunique(), df_encode['tf'].nunique(), np.union1d(df_jaspar['tf'].unique(), df_encode['tf'].unique()).shape) #634 957 (1282,)
            print(df_jaspar['peaks'].nunique(), df_encode['peaks'].nunique(), np.union1d(df_jaspar['peaks'].unique(), df_encode['peaks'].unique()).shape) #62310613 19645880 (81956493,)
        df_concat = pd.concat([df_jaspar, df_encode], axis=0)
        df_concat = df_concat[cols].drop_duplicates() #(143124468, 5)
        df_concat.to_csv('output/db/jaspar_encode.bed.gz', sep='\t', header=False, index=False, compression='gzip')


def process_trace_seqera(trace):
    trace['model'] = pd.DataFrame(trace.name.str.split(':').to_list())[3] #TODO: number 3 might be different  
    trace = trace.groupby('model').apply(lambda df: df.sort_values(by='duration', ascending=False).iloc[0])[['%cpu', 'peak_rss', 'peak_vmem', 'rchar', 'wchar', 'duration']]

    def convert_duration_to_hours(duration_str):
        import re
        hours, minutes, seconds = 0, 0, 0
        time_parts = re.findall(r'(\d+)([hms])', duration_str)
        for value, unit in time_parts:
            if unit == 'h':
                hours = int(value)
            elif unit == 'm':
                minutes = int(value)
            elif unit == 's':
                seconds = int(value)
        return (hours * 3600 + minutes * 60 + seconds)/3600
    def format_ram(row):
        value = float(row.split()[0])
        if 'GB' in row:
            value = value
        elif 'MB' in row:
            value = value/1000
        else:
            raise ValueError('Define')
        return value 

    for col in trace.columns:
        if col=='%cpu':
            trace[col] = trace[col].str.replace('%', '').astype(float)
        elif col=='duration':
            trace[col] = trace[col].apply(convert_duration_to_hours)
        else:
            trace[col] = trace[col].apply(format_ram)
    return trace 



if __name__ == '__main__':
    # run_grn_inference()
    # calculate_scores()
    # analyse_meta_cells(task_grn_inference_dir='./')
    analyse_imputation(task_grn_inference_dir='./')
    # analyse_corr_vs_tfmasked_corr(task_grn_inference_dir='./')

