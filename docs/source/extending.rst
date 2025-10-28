Extending
======================
In this section, we provide instructions on how to add new GRN inference methods, metrics, and datasets to the geneRNIB platform. 
Installing geneRNIB is a prerequisite for integrating new features. 

Add a GRN inference method
-----------------------------------

Examples of GRN inference methods include GRNBoost2, CellOracle, and SCENIC. 
The list of integrated GRN inference methods can be found on the geneRNIB platform, `src/methods`, which are examples of how to integrate new methods for both R and Python. 

Each method requires a `config.vsh` file together
with a `script.py`. Additionally, the method can have extra files to store and organize the code, such as `helper`, which are stored in the same folder and called by the `script.py`.

The overlook of `config.vsh` is as follows. However, refer to the `src/methods/` folder for the updated formatting.

.. code-block:: yaml
   :caption: Example of a `config.vsh` file

   __merge__: /src/api/comp_method.yaml # merge with the common method scheme

    name: method_name # unique id for your method
    namespace: "grn_methods"
    info:
    label: pretty method name # e.g. "GRNBoost2"
    summary: "summary of your method"
    description: |
       more about your method 
    documentation_url: link to your method documentation # optional

    resources:
        - type: python_script # or R_script
          path: script.py # your main script (dont change it). or script.R

    engines:
        - type: docker 
            image: ghcr.io/openproblems-bio/base_python:1.0.4 # a base docker image
            __merge__: /src/api/base_requirements.yaml # merge with the base requirements schema required for the pipeline
            setup:
            - type: python
                packages: [ grnboost2 ] # additional packages required for your method. see different methods for examples as this could get complicated. or, use your image and omit this.
    
    runners: # this is for the nextflow pipeline.
        - type: executable
        - type: nextflow
            directives:
            label: [midtime, midmem, midcpu] # expected resources. see scripts/labels_tw.config for their definition 

Your `script.py` should have the following structure:

.. code-block:: python
   :caption: Example of a `script.py` file

   import sys 
   import anndata as ad
   import numpy as np
   ... # import necessary libraries

    ## VIASH START  -> this is necessary for the viash to work. It essentially replaces this with the parameters passed to the config.vsh file
    par = {
        'rna': 'resources/grn_benchmark/rna_op.h5ad', # example rna data
        'tf_all': 'resources/grn_benchmark/prior/tf_all.csv', # tf list. you will need this to filter the network to only tf-gene pairs. we only evaluate top 50k TF-gene edges so better to filter it.
        'prediction': 'output/prediction.h5ad' # where to save the prediction
    }
    ## VIASH END

    # Load some data
    rna = ad.read_h5ad(par["rna"])
    tf_all = np.loadtxt(par["tf_all"], dtype=str)

    # Your method code here
    net = pd.DataFrame({"source": ["TF1", "TF2"], "target": ["gene1", "gene2"], "weight": [0.1, 0.2]}) # example network


    # Save the inferred network
    net['weight'] = net['weight'].astype(str)  # Ensure weight is stored as a string
    output = AnnData(
        X=None,
        uns={
            "method_id": "method_name",
            "dataset_id": "dataset_name", 
            "prediction": net[["source", "target", "weight"]]
        }
    )
    output.write(params["prediction"])

Once you have added your method, you can test it by running the following command. For this, download and place the test datasets in `resources_test/grn_benchmark`.

.. code-block:: bash
    :caption: Download test data

    aws s3 sync s3://openproblems-data/resources_test/grn/grn_benchmark resources_test/grn_benchmark --no-sign-request 


.. code-block:: bash
   :caption: Test your method

   viash test src/methods/your_method/config.vsh # path to the config.vsh file of your method

Once the test is successful, you can submit a pull request to the geneRNIB repository to integrate your method.
See additional Viash commands in the `Viash documentation <https://viash.io/guide/>`_ to run your method with different parameters.

Add a GRN evaluation metric
-----------------------------------
Similar to method integration, metrics also follow similar file formatting. See folder `src/metrics/` for examples.
While new metrics could use different evaluation datasets, the current files for evaluation are located in `resources/grn_benchmark/evaluation_data`.
There are three formats of datasets; single cell, (pseudo)bulk, and differential expression (de). 
The choice of the dataset depends on the evaluation metric.

A few tips:

- use `read_prediction` from from `src/utils/util` to read the inferred GRNs and do checking the format.
- the metric should output a score which has to be a h5ad file. In Python, this can be done:

.. code-block:: python

    import pandas as pd
    from util import format_save_score

    results = pd.DataFrame({
        'metric_key_1': [metric_value_1],  # submetric 1
        'metric_key_2': [metric_value_2],  # submetric 2
    })
    method_id = 'name of GRN method'
    dataset_id = 'name of dataset used for GRN inference'
    score_file = 'output/score.h5ad'

    format_save_score(results, method_id, dataset_id, score_file)

Add a GRN inference and evalaution dataset
-----------------------------------
Here we explain how to integrate new datasets. All datasets are in h5ad, and the example structure of a inference or evaluation dataset can be found in `resources/grn_benchmark/`.
The inference datasets are in `resources/grn_benchmark/inference_data/` and the evaluation datasets are in `resources/grn_benchmark/evaluation_data/`.
Each dataset should have a unique `dataset_id`, stored in `.uns['dataset_id']`, that will be used to identify it in the platform.
In addition, there should be additional information stored in `.uns`, such as the description, reference, and normalization type.
Also, the normalized values should be stores in `layers` with the name of normalization method, e.g. lognorm.
