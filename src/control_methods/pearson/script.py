
## VIASH START
par = {
    'multiomics_rna': 'resources/grn-benchmark/multiomics_rna_0.h5ad',
    'tf_all': 'resources/prior/tf_all.csv',
    'cell_type_specific': True,
    'max_n_links': 50000,
    'prediction': 'resources/grn_models/donor_0_default/pearson.csv',
    "seed": 32,
    'only_hvgs': True,
    'normalize': False
}
## VIASH END
import sys
sys.path.append(meta["resources_dir"])
from util import create_corr_net

par['causal'] = False
net = create_corr_net(par)
print('Output GRN')
net.to_csv(par['prediction'])
