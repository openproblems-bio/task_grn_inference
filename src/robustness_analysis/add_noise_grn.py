import os
import pandas as pd
import numpy as np

layer = 'scgen_pearson'
grn_folder = 'resources/grn_models'
grn_folder_noised = 'resources/supplementary/grn_models_noised'
noise_ratio = 0.2

# Ensure the output folder exists
os.makedirs(grn_folder_noised, exist_ok=True)

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
