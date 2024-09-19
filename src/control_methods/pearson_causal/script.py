
## VIASH START
par = {
    'multiomics_rna': 'resources/grn-benchmark/multiomics_rna_0.h5ad',
    'tf_all': 'resources/prior/tf_all.csv',
    'cell_type_specific': True,
    'max_n_links': 50000,
    'prediction': 'resources/grn_models/donor_0_default/pearson_causal.csv',
    "seed": 32
}
## VIASH END
# import sys
# sys.path.append('./src/utils')
from util import create_corr_net

print('Create causal corr net')
par['causal'] = True
net = create_corr_net(par)

print('Output GRN')
net.to_csv(par['prediction'])
