Run GRN Inference
=================

GRN Inference Without Method Integration
----------------------------------------

In this section, we explain how to access datasets and infer gene regulatory networks (GRNs) using your method without integrating it into geneRNIB.

### 1. Download the Inference Datasets  
The inference datasets are stored on AWS and can be downloaded using the following command:

.. code-block:: bash

   aws s3 sync s3://openproblems-data/resources/grn/grn_benchmark/inference_data resources/grn_benchmark/inference_data --no-sign-request

### 2. Available Datasets  
The available datasets include **op, nakatake, replogle, adamson,** and **norman**. Each dataset provides RNA data. The `op` dataset also includes ATAC data.

### 3. GRN Inference Guidelines  
When performing GRN inference, consider the following:  

- We evaluate only the **top TF-gene pairs**, currently limited to **50,000 edges**, ranked by their assigned weight.  
- The inferred network should follow this format:  

  **Columns:**  
  - `source`: Transcription factor (TF)  
  - `target`: Target gene  
  - `weight`: Regulatory importance/likelihood score/etc.  

### 4. Saving the Inferred Network  
Since geneRNIB works with **AnnData**, your inferred network should be saved in this format.

If your network is a pandas DataFrame with three columns (`source`, `target`, `weight`), you can save it as follows (replace `grnboost2` and `norman` with your method and the dataset name used to infer the GRN):

.. code-block:: python
   net['weight'] = net['weight'].astype(str)  # Ensure weight is stored as a string

   output = ad.AnnData(
       X=None,
       uns={
           "method_id": "grnboost2",
           "dataset_id": "norman",
           "prediction": net[["source", "target", "weight"]]
       }
   )

   output.write("save_to_file.h5ad")

For R, use the following approach:
.. 
.. code-block:: r
   net$weight <- as.character(net$weight)

   output <- AnnData(
       X = matrix(nrow = 0, ncol = 0),
       uns = list(
           method_id = "grnboost2",
           dataset_id = "norman",
           prediction = net[, c("source", "target", "weight")]
       )
   )

   output$write_h5ad("save_to_file.h5ad", compression = "gzip")

### Next Steps  
Once you have inferred GRNs for one or more datasets, proceed to the next section to run the evaluation.
