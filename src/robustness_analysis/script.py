import os
import pandas as pd
import numpy as np

## VIASH START
par = {
  "prediction": "resources/grn_models/collectri.csv",
  "prediction_n": "output/grn_noised.csv",
  'degree': 20,
  'type': 'links'
}

## VIASH END

degree = par['degree']/100

prediction = pd.read_csv(par['prediction'])
assert 'weight' in prediction.columns 

if type =='weight':
    print('Add noise to weight')
    std_dev = prediction['weight'].std()
    noise = np.random.normal(0, degree * std_dev, size=prediction['weight'].shape)
    prediction['weight'] += noise

elif type =='links':
    print('Permute links')
    num_rows_to_permute = int(len(prediction) * degree)
    permute_indices = np.random.choice(prediction.index, size=num_rows_to_permute, replace=False)
    
    prediction.loc[permute_indices, 'weight'] = np.random.permutation(prediction.loc[permute_indices, 'weight'].values)
print('Output noised GRN')
prediction.to_csv(par['prediction_n'])

