import os
import pandas as pd
import numpy as np

layer = 'scgen_pearson'
grn_folder = 'resources/grn_models'
grn_folder_noised = 'resources/supplementary/grn_models_noised'
noise_ratio = 0.2
# permute_ratio = 0.2

# Ensure the output folder exists
os.makedirs(grn_folder_noised, exist_ok=True)

if True: # add noise
    # Loop through all files in the grn_folder
    for file_name in os.listdir(grn_folder):
        if file_name.endswith('.csv'):
            # Read the CSV file
            file_path = os.path.join(grn_folder, file_name)
            df = pd.read_csv(file_path)

            # Add noise to the 'weight' column
            if 'weight' in df.columns:
                std_dev = df['weight'].std()
                noise = np.random.normal(0, noise_ratio * std_dev, size=df['weight'].shape)
                df['weight'] += noise

            # Save the noised DataFrame to the new folder
            noised_file_path = os.path.join(grn_folder_noised, file_name)
            df.to_csv(noised_file_path, index=False)

    print("Noise added to all GRN models and saved successfully.")
# Loop through all files in the grn_folder
else:
    for file_name in os.listdir(grn_folder):
        if file_name.endswith('.csv'):
            # Read the CSV file
            file_path = os.path.join(grn_folder, file_name)
            df = pd.read_csv(file_path)

            # Permute 20% of the rows in the 'weight' column
            if 'weight' in df.columns:
                num_rows_to_permute = int(len(df) * permute_ratio)
                
                # Randomly select 20% of the row indices to permute
                permute_indices = np.random.choice(df.index, size=num_rows_to_permute, replace=False)
                
                # Shuffle the selected rows in 'weight' column
                df.loc[permute_indices, 'weight'] = np.random.permutation(df.loc[permute_indices, 'weight'].values)

            # Save the modified DataFrame to the new folder
            noised_file_path = os.path.join(grn_folder_noised, file_name)
            df.to_csv(noised_file_path, index=False)

    print("20% of the 'weight' column rows have been permuted for all GRN models and saved successfully.")