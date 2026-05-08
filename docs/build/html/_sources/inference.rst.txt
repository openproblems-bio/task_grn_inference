GRN Inference
=================

GRN inference using integrated methods
----------------------------------------

To infer a GRN for a given dataset using an integrated method, first download the test datasets:

.. code-block:: bash

   aws s3 sync s3://openproblems-data/resources_test/grn/grn_benchmark \
     resources_test/grn_benchmark --no-sign-request

Use ``resources_test/`` paths for development and testing. Switch to ``resources/grn_benchmark/`` for full-scale runs (requires downloading the main datasets — see :doc:`dataset`).

**With Docker (Viash)**

.. code-block:: bash

   viash run src/methods/pearson_corr/config.vsh.yaml -- \
     --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
     --prediction output/net.h5ad \
     --tf_all resources_test/grn_benchmark/prior/tf_all.csv

**Without Docker (conda or Singularity)**

.. code-block:: bash

   bash src/methods/pearson_corr/run_local.sh \
     --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
     --tf_all resources_test/grn_benchmark/prior/tf_all.csv \
     --prediction output/net.h5ad

Replace ``pearson_corr`` with any method under ``src/methods/``. Methods that require Singularity (e.g. GRNBoost2, scPRINT, Scenic/Scenic+, CellOracle) will invoke it automatically from their ``run_local.sh`` script.


GRN inference without method integration
----------------------------------------

Download the inference datasets from the :doc:`dataset` page. For each dataset, perform GRN inference using your own algorithm. Consider the following:

- We evaluate only the **top TF-gene pairs**, currently limited to **50,000 edges**, ranked by their assigned weight.
- The inferred network should follow this format:

  **Columns:**

  - ``source``: Transcription factor (TF)
  - ``target``: Target gene
  - ``weight``: Regulatory importance/likelihood score

Since geneRNIB works with **AnnData**, your inferred network should be saved in this format.
If your network is a pandas DataFrame with three columns (``source``, ``target``, ``weight``), you can save it as follows:

**Python**

.. code-block:: python

   import anndata as ad

   net['weight'] = net['weight'].astype(str)  # weight must be stored as a string

   output = ad.AnnData(
       X=None,
       uns={
           "method_id": "your_method_name",  # e.g. "grnboost2"
           "dataset_id": "dataset_name",     # e.g. "replogle"
           "prediction": net[["source", "target", "weight"]]
       }
   )
   output.write_h5ad("save_to_file.h5ad")

**R**

.. code-block:: r

   library(zellkonverter)

   net$weight <- as.character(net$weight)

   output <- AnnData(
       X = matrix(nrow = 0, ncol = 0),
       uns = list(
           method_id = "your_method_name",  # e.g. "grnboost2"
           dataset_id = "dataset_name",     # e.g. "replogle"
           prediction = net[, c("source", "target", "weight")]
       )
   )
   output$write_h5ad("save_to_file.h5ad", compression = "gzip")

Once you have inferred GRNs for one or more datasets, proceed to the :doc:`evaluation` page to run the evaluation metrics.
