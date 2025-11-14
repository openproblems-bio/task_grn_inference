GRN Inference
=================

GRN inference using integrated methods
----------------------------------------
To infer a GRN for a given dataset (e.g. op) using our integration methods, after installation of geneRNBI, run the following command (e.g. Pearson correlation):

.. code-block:: bash

  viash run src/methods/pearson_corr/config.vsh.yaml -- \
    --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
    --prediction output/net.h5ad \
    --tf_all resources_test/grn_benchmark/prior/tf_all.csv

See `src/methods/` for the list of currently integrated methods.


GRN inference without method integration
----------------------------------------
Download the inference datasets from the `dataset` page. For each dataset, perform GRN inference using your algorithm. Consider the following:  

- We evaluate only the **top TF-gene pairs**, currently limited to **50,000 edges**, ranked by their assigned weight.  
- The inferred network should follow this format:  

  **Columns:**  
  
  - `source`: Transcription factor (TF)  
  - `target`: Target gene  
  - `weight`: Regulatory importance/likelihood score  

Since geneRNIB works with **AnnData**, your inferred network should be saved in this format.
See `resources/grn_benchmark/prior/collectri.h5ad` for an example of the expected format.
If your network is a pandas DataFrame with three columns (`source`, `target`, `weight`), you can save it as follows:

.. code-block:: python

   import anndata as ad

   net['weight'] = net['weight'].astype(str)  # Ensure weight is stored as a string

   output = ad.AnnData(
       X=None,
       uns={
           "method_id": "your_method_name", # e.g. "grnboost2"
           "dataset_id": "dataset_name",  #e.g. "replogle"
           "prediction": net[["source", "target", "weight"]]
       }
   )

   output.write("save_to_file.h5ad")

For R, use the following approach:

.. code-block:: r

   library(zellkonverter)

   net$weight <- as.character(net$weight)

   output <- AnnData(
       X = matrix(nrow = 0, ncol = 0),
       uns = list(
           method_id = "your_method_name", # e.g. "grnboost2"
           dataset_id ="dataset_name",  # e.g. "replogle"
           prediction = net[, c("source", "target", "weight")]
       )
   )

   output$write_h5ad("save_to_file.h5ad", compression = "gzip")

Once you have inferred GRNs for one or more datasets, proceed to the next section to run the evaluation.
