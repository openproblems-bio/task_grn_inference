import pandas as pd
import numpy as np
import mudata as mu
import os
import argparse
import warnings
warnings.filterwarnings("ignore")


parser = argparse.ArgumentParser(description="", usage="")
parser.add_argument('--pre_path', required=True)
parser.add_argument('--exp_path', required=True)
parser.add_argument('--pks_path', required=True)

args = vars(parser.parse_args())
pre_path = args['pre_path']
exp_path = args['exp_path']
pks_path = args['pks_path']


# Write the RNA matrix
pre_type = os.path.basename(pre_path).split('.')[0]
data = mu.read(pre_path)
print(data['rna'].obs)
rna_X = pd.DataFrame(np.array(data['rna'].layers['counts'].todense()).T, columns=data['rna'].obs.index, index=data['rna'].var.index)
rna_X.to_csv(exp_path, sep="\t", compression="gzip")

all_atac_peak = np.unique([n.replace(':', '-') for n in data['atac'].var.index])

all_atac_peak = pd.DataFrame([n.split('-') for n in all_atac_peak])
all_atac_peak.columns = ['chr', 'srt', 'end']
all_atac_peak['srt'] = all_atac_peak['srt'].astype(int)
all_atac_peak['end'] = all_atac_peak['end'].astype(int)
all_atac_peak = all_atac_peak[(all_atac_peak.end - all_atac_peak.srt) >= 100]
all_atac_peak = all_atac_peak.sort_values(by=['chr', 'srt', 'end'])
all_atac_peak.to_csv(pks_path, sep='\t', header=False, index=False)

# Store clusters
clus = sorted(data['rna'].obs['cell_type'].unique())
for c in clus:
    ctype_ids = data[data['rna'].obs['cell_type'] == c].obs.index
    c = c.replace(' ', '_')
    with open(os.path.join(os.path.dirname(exp_path), f'barcodes_{c}.txt'), "w") as f:
        for i in ctype_ids:
            f.write(f"{i}\n")
