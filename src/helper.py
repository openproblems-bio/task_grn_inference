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

def run_grn_inference():
    # par = {
    #     'methods': ['portia'],
    #     'models_dir': 'resources/grn_models/mccalla/han',
    #     'multiomics_rna': 'resources/grn-benchmark/mccalla/inference/han.h5ad', 
    #     'num_workers': 20,
    #     'mem': "250GB",
    #     'time': "48:00:00",
    #     'max_n_links': 100000,
    #     'causal': False, 
    #     'normalize': False
    #     }
    par = {
        'methods': ['scenicplus'],
        'models_dir': 'resources/grn_models/',
        'multiomics_rna': 'resources/grn-benchmark/multiomics_rna.h5ad', 
        'multiomics_atac': 'resources/grn-benchmark/multiomics_atac.h5ad', 
        'num_workers': 20,
        'mem': "400GB",
        'time': "48:00:00",
        'causal': False,
        'normalize': True,
        'max_n_links': 100000,
    }

    for method in par['methods']:
        print(method)
        par['prediction'] = f"{par['models_dir']}/{method}.csv"
        
        # Method arguments
        method_args = (f"--multiomics_rna {par['multiomics_rna']} "
                       f"--prediction {par['prediction']} "
                       f"--num_workers {par['num_workers']} "
                       f"--max_n_links {par['max_n_links']} ")
        if par['causal']:
            method_args += f"--causal "

        # Determine the command based on the method
        if method in ["pearson_corr", "positive_control", "negative_control"]:
            if par['normalize']:
                method_args += f"--normalize "
            command = f"python src/control_methods/{method}/script.py {method_args}"
        elif method == "celloracle":
            method_args += f"--multiomics_atac {par['multiomics_atac']} "
            command = (f"/home/jnourisa/miniconda3/envs/celloracle/bin/python "
                       f"src/methods/multi_omics/celloracle/script.py {method_args}")
        elif method in ["grnboost2", "scenic", "genie3"]:
            command = f"singularity exec ../../images/scenic python src/methods/single_omics/{method}/script.py {method_args}"
        elif method == 'scglue':
            method_args += f"--multiomics_atac {par['multiomics_atac']} "
            command = f"singularity exec ../../images/scglue python src/methods/multi_omics/{method}/script.py {method_args}"
        elif method == 'scenicplus':
            method_args += f"--multiomics_atac {par['multiomics_atac']} "
            command = f"singularity exec ../../images/scenicplus python src/methods/multi_omics/{method}/script.py {method_args}"
        else:
            command = f"singularity exec ../../images/{method} python src/methods/single_omics/{method}/script.py {method_args}"

        # Prepare sbatch command
        tag = f"--job-name={method}"  # No spaces around '='
        resources = (f"--cpus-per-task={par['num_workers']} "
                     f"--mem={par['mem']} --time={par['time']}")
        
        # Combine tags and resources
        full_tag = [tag] + resources.split()
        
        # Add GPU partition if method is 'scglue'
        if method == 'scglue':
            full_tag += ["--partition=gpu", "--gres=gpu:1"]
        if True:
            full_tag += ["--partition=gpu", "--gres=gpu:1"]
        # Run sbatch command
        try:
            result = subprocess.run(['sbatch'] + full_tag + ['scripts/sbatch/grn_inference.sh', command], check=True, capture_output=True, text=True)
            # result = subprocess.run(['bash'] + ['scripts/sbatch/grn_inference.sh', command], check=True, capture_output=True, text=True)

            print(f"Job {method} submitted successfully.")
            print(result.stdout)  # Print the standard output
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while submitting job for {method}: {e}")
            print(f"Return code: {e.returncode}")
            print(f"Command output: {e.output}")
            print(f"Command error output: {e.stderr}")


def create_skeleton():
    # - run src/metrics/skeleton for each motif 
    # - combined the resulting tf2gene
    
    skeleton_encode = pd.read_csv('output/skeleton_encode/tf2gene.csv', index_col=0)
    skeleton_jaspar = pd.read_csv('output/skeleton_jaspar/tf2gene.csv', index_col=0)
    skeleton = pd.concat([skeleton_encode, skeleton_jaspar], axis=0)
    skeleton = skeleton.drop_duplicates().reset_index(drop=True)
    skeleton.nunique()
    all_links = skeleton['source'].astype(str) + '_' + skeleton['target'].astype(str)
    np.savetxt('resources/prior/skeleton.csv', all_links.values, fmt='%s')

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

def marco_data():
    # import subprocess
    # import anndata as ad 
    # import pandas as pd
    # import numpy as np
    # for cell_type in ['zhao', 'shalek', 'han', 'jackson']:
    #     adata = ad.read_h5ad(f'resources_local/mccalla_extended/{cell_type}.h5ad')
    #     adata.layers['norm'] = adata.X
    #     adata.obs['cell_type'] = 'onecelltype'
    #     adata.write(f'resources_local/mccalla_extended/{cell_type}.h5ad')
    #     subsample = min([10000, len(adata)])
    #     for GT in ['KDunion', 'chipunion', 'chipunion_KDUnion_intersect']:
    #         GT_df = pd.read_csv(f'resources_local/mccalla_extended/{cell_type}_{GT}.csv')
    #         gene_overlap = np.intersect1d(adata.var_names, GT_df.target.unique()).shape
    #         print(f"{cell_type}-{GT}. adata shape: {adata.shape}, GT size: {GT_df.shape}, Gene overlap: {gene_overlap}")
    #         command = f"viash run src/metrics/regression_1/config.vsh.yaml -- --perturbation_data resources_local/mccalla_extended/{cell_type}.h5ad --prediction resources_local/mccalla_extended/{cell_type}_{GT}.csv --layer norm --subsample {subsample} --apply_tf false --tf_all resources/prior/tf_all.csv --max_n_links -1 --verbose 1 --score output/{cell_type}_{GT}.h5ad"
    #         subprocess.run(command, shell=True, check=True)


    all_gene_names = pd.read_csv('resources/prior/genome_annotation.tsv', sep='\t').Gene.unique().astype(str)

    # adata_names = ['han', 'jackson', 'zhao', 'shalek']
    # post_fixes = ['chipunion','KDunion', 'chipunion_KDUnion_intersect']
    # for name in adata_names:
    #     print('------- ', name)
    #     adata = ad.read_h5ad(f'resources/grn-benchmark/mccalla/inference/{name}.h5ad')
    #     adata.var_names = adata.var_names.str.upper().astype(str)

    #     genes = adata.var_names
    #     print(len(genes),genes.isin(all_gene_names).sum())
    #     print(np.setdiff1d(genes, all_gene_names)[0:5])
    #     adata.write_h5ad(f'resources/grn-benchmark/mccalla/inference/{name}.h5ad')

    #     for post_fix in post_fixes:
    #         GT = pd.read_csv(f'resources/grn-benchmark/mccalla/evaluation/{name}_{post_fix}.csv')
    #         GT.source = GT.source.str.upper()
    #         GT.target = GT.target.str.upper()
    #         GT.to_csv(f'resources/grn-benchmark/mccalla/evaluation/{name}_{post_fix}.csv')
    #         tf_genes =set(GT.source) | set(GT.target)
    #         tf_genes = [name.upper() for name in tf_genes]
    #         # print(tf_genes)
    #         print('-- ', post_fix, len(tf_genes), np.intersect1d(list(tf_genes), genes).shape)
    pass

def extract_data(data, reg='reg1', dataset_id='scgen_pearson'):
    i = 0
    for entry in data:
        if entry['dataset_id']!=dataset_id:
            continue
        try:
            rg, method_id = entry['method_id'].split('-')
        except:
            rg, method_id, _ = entry['method_id'].split('-')
        if rg != reg:
            continue
        dataset_id = entry['dataset_id']
        metric_ids = entry['metric_ids']
        metric_values = entry['metric_values']
        
        df = pd.DataFrame([metric_values], index=[method_id], columns=metric_ids)
        if i==0:
            df_reg = df
        else:
            df_reg = pd.concat([df_reg, df], axis=0)
        i+=1
    return df_reg
def process_data(RUN_ID, models_all=None):
    raise ValueError('fix this')
    # !aws s3 sync s3://openproblems-data/resources/grn/results/{RUN_ID} resources/results/{RUN_ID} 
    base_folder = f'resources/results/{RUN_ID}'
    result_file = f'{base_folder}/scores.yaml'
        

    with open(result_file, 'r') as file:
        data = yaml.safe_load(file)
    
    
    if models_all is None:
        df_reg1 = extract_data(data, reg='reg1')
        df_reg2 = extract_data(data, reg='reg2')
    else:
        df_reg1 = extract_data(data, reg='reg1').reindex(models_all)
        df_reg2 = extract_data(data, reg='reg2').reindex(models_all)
    # df_all = pd.concat([df_reg1, df_reg2], axis=1).fillna(0)
    # df_all_n = (df_all-df_all.min(axis=0))/(df_all.max(axis=0)-df_all.min(axis=0))
    # df_all['Rank'] = df_all_n.mean(axis=1).rank(ascending=False).astype(int)
    df_all = pd.concat([df_reg1, df_reg2], axis=1)
    return df_all
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
def process_trace_local(job_ids_dict):
    def get_sacct_data(job_id):
        command = f'sacct -j {job_id} --format=JobID,JobName,AllocCPUS,Elapsed,State,MaxRSS,MaxVMSize'
        output = subprocess.check_output(command, shell=True).decode()
        
        # Load the output into a DataFrame
        df = pd.read_csv(io.StringIO(output), delim_whitespace=True)
        df = df.iloc[[2]]
        return df
    def elapsed_to_hours(elapsed_str):
        time = elapsed_str.split('-')
        if len(time) > 1:
            day = int(time[0])
            time = time[1]
        else:
            day = 0
            time = time[0]
        h, m, s = map(int, time.split(':'))
        return day*24 + h + m / 60 + s / 3600
    def reformat_data(df_local):
        # Remove 'K' and convert to integers
        df_local['MaxRSS'] = df_local['MaxRSS'].str.replace('K', '').astype(int)
        df_local['MaxVMSize'] = df_local['MaxVMSize'].str.replace('K', '').astype(int)
        df_local['Elapsed'] = df_local['Elapsed'].apply(lambda x: (elapsed_to_hours(x)))

        # Convert MaxRSS and MaxVMSize from KB to GB
        df_local['MaxRSS'] = df_local['MaxRSS'] / (1024 ** 2)  # Convert KB to GB
        df_local['MaxVMSize'] = df_local['MaxVMSize'] / (1024 ** 2)  # Convert KB to GB
        return df_local
    for i, (name, job_id) in enumerate(job_ids_dict.items()):
        if type(job_id)==list:
            
            for i_sub, job_id_ in enumerate(job_id):
                df_ = get_sacct_data(job_id_)
                df_ = reformat_data(df_)
                if i_sub == 0:
                    df = df_
                else:
                    concat_df = pd.concat([df, df_], axis=0)
                    df['MaxVMSize'] = concat_df['MaxVMSize'].max()
                    df['MaxRSS'] = concat_df['MaxRSS'].max()
                    df['Elapsed'] = concat_df['Elapsed'].sum()
        else: 
            df = get_sacct_data(job_id)
            df = reformat_data(df)
        df.index = [name]
        if i==0:
            df_local = df
        else:
            df_local = pd.concat([df_local, df], axis=0)
        
    
    return df_local


if __name__ == '__main__':
    run_grn_inference()
    # calculate_scores()