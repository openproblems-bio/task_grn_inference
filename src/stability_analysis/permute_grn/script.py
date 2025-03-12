import os
import pandas as pd
import numpy as np
import sys
import anndata as ad


## VIASH START
par = {
  "prediction": "resources/grn_models/collectri.h5ad",
  "prediction_n": "output/grn_noised.h5ad",
  "noise_type": "weight",
  "degree": 10
}
## VIASH END


def main(par):
  net = ad.read_h5ad(par['prediction'])
  prediction = pd.DataFrame(net.uns['prediction'])
  prediction['weight'] = prediction['weight'].astype(float)

  degree = par['degree']/100
  type = par['noise_type']

  if type == 'weight': # add noise to weight
    assert 'weight' in prediction.columns 
    print('Add noise to weight')
    std_dev = prediction['weight'].std()
    noise = np.random.normal(loc=0, scale=degree * std_dev, size=prediction['weight'].shape)
    prediction['weight'] += noise

  elif type == 'net': # shuffle source-target matrix
    print('Permute links')
      
    # 1. Pivot the GRN with target as index and source as columns
    prediction = prediction.groupby(['target', 'source'], as_index=False)['weight'].mean()

    pivot_df = prediction.pivot(index='target', columns='source', values='weight')

    # Fill NaNs with 0 (no interaction)
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
  elif type == 'direction': # change the regulatory sign
    # Calculate the number of rows to permute
    prediction = prediction.reset_index(drop=True)
    n_rows_to_permute = int(len(prediction) * (degree))
    # print(n_rows_to_permute)
    
    # Randomly select indices to permute
    indices_to_permute = np.random.choice(prediction.index, size=n_rows_to_permute, replace=False)

    print(indices_to_permute)
    # Swap source and target for the selected rows
    prediction.loc[indices_to_permute, ['source', 'target']] = prediction.loc[indices_to_permute, ['target', 'source']].values
    
  else:
    raise ValueError(f'Wrong type ({type}) for adding noise')
  
  prediction = prediction.astype(str)
  net = ad.AnnData(X=None, uns={"method_id": net.uns['method_id'], "dataset_id": net.uns['dataset_id'], "prediction": prediction[["source", "target", "weight"]]})
  net.write(par['prediction_n'])

if __name__ == '__main__':
  main(par)

