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
print('Output noised GRN')
prediction = main(par)
prediction.to_csv(par['prediction_n'])

