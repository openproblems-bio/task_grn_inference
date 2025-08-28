import yaml
import os
import pandas as pd
import anndata as ad 
import numpy as np
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



# if __name__ == '__main__':
    # run_grn_inference()
    # calculate_scores()
    # analyse_meta_cells(task_grn_inference_dir='./')
    # analyse_imputation(task_grn_inference_dir='./')
    # analyse_corr_vs_tfmasked_corr(task_grn_inference_dir='./')

