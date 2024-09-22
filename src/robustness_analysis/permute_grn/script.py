import os
import pandas as pd
import numpy as np

## VIASH START
par = {
  "prediction": "resources/grn_models/collectri.csv",
  "prediction_n": "output/grn_noised.csv",
  'degree': 20,
  'noise_type': 'links'
}

## VIASH END

degree = par['degree']/100
type = par['noise_type']


prediction = pd.read_csv(par['prediction'])


if type == 'weight': # add noise to weight
  assert 'weight' in prediction.columns 
  print('Add noise to weight')
  std_dev = prediction['weight'].std()
  noise = np.random.normal(loc=0, scale=degree * std_dev, size=prediction['weight'].shape)
  prediction['weight'] += noise

elif type == 'net': # shuffle source-target matrix
  print('Permute links')
    
  # 1. Pivot the GRN with target as index and source as columns
  pivot_df = prediction.pivot(index='target', columns='source', values='weight')

  # Fill NaNs with 0 or a value of your choice
  pivot_df.fillna(0, inplace=True)

  # 2. Randomly choose degree% of the matrix to shuffle
  matrix_flattened = pivot_df.values.flatten()
  n_elements = len(matrix_flattened)
  n_shuffle = int(n_elements * degree)

  # Randomly select 20% of the matrix elements' indices
  shuffle_indices = np.random.choice(n_elements, n_shuffle, replace=False)

  # Get the values that will be shuffled
  shuffle_values = matrix_flattened[shuffle_indices]

  # 3. Shuffle the selected values
  np.random.shuffle(shuffle_values)

  # Assign the shuffled values back to the selected positions
  matrix_flattened[shuffle_indices] = shuffle_values

  # Reshape the flattened array back into the matrix
  pivot_df_shuffled = pd.DataFrame(matrix_flattened.reshape(pivot_df.shape), 
                                  index=pivot_df.index, 
                                  columns=pivot_df.columns)
                              
  flat_df = pivot_df_shuffled.reset_index()

  # Melt the DataFrame to turn it back into long-form (source-target-weight)
  prediction = flat_df.melt(id_vars='target', var_name='source', value_name='weight')
  prediction = prediction[prediction['weight'] !=0 ].reset_index(drop=True)
elif type == 'sign': # change the regulatory sign
  num_rows = len(prediction)
  num_to_modify = int(num_rows * degree)
  # 2. Randomly select indices to modify
  random_indices = np.random.choice(prediction.index, size=num_to_modify, replace=False)
  # 3. Change the sign of the selected rows
  prediction.loc[random_indices, 'weight'] *= -1
elif type == 'binary': # change the regulatory sign
  prediction['weight'] = np.where(prediction['weight'] > 0, 1, -1)
else:
  raise ValueError(f'Wrong type ({type}) for adding noise')

print('Output noised GRN')
prediction.to_csv(par['prediction_n'])

