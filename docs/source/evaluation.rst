
GRN evaluation
=================
The evaluation metrics used in geneRNIB are summarized below. For a detailed description of each metric, refer to the geneRNIB paper.

We originally defined **eight evaluation metrics**, grouped into three categories: **Regression 1, Regression 2, and Wasserstein Distance**. 
However, we recently removed **Regression 1** as it did not prove to be effective for perturbational settings. 

- The **regression-based metrics** assess the predictive power of an inferred GRN by using regression models to predict perturbation data (evaluation data) based on the feature space constructed from the inferred network.  
- The **Wasserstein distance-based metric** evaluates GRN edges by measuring the distributional shift in target gene expression between observations and perturbation data for a given transcription factor (TF).

Wasserstein distance-based metrics are only applicable for datasets that are gene perturbations and are in single cell format. Thus, currently the following datasets are supported:
- Replogle
- Xaira:HEK293T
- Xaira:HCT116
- Norman
- Adamson
  
.. image:: images/metrics.png
   :width: 90%
   :align: center
----

The evaluation metrics expect the inferred network to be in the form of an AnnData object with specific format as explained here. It should be noted that the metric currently evaluate only the **top TF-gene pairs**, currently limited to **50,000 edges**, ranked by their assigned weight.  

The inferred network should have a tabular format with the following columns:  

  - `source`: TF gene name
  - `target`: Target gene gene  
  - `weight`: Regulatory importance/likelihood score/etc.  

See `resources_test/grn_models/op/collectri.h5ad` for an example of the expected format.

For the regression based approaches, we used the pseudobulk version of the perturbation data while for the Wasserstein distance, the single cell data are used.

It should be noted that for Wasserstein distance, we have already computed all possible combination of TF-gene pairs and stored it in the `resources/grn_benchmark/prior/` folder.
This substantially reduces the computation time during evaluation.

To run the evalution for a given GRN and dataset, use the following command:
```bash
bash scripts/run_grn_evaluation.sh --prediction=<inferred GRN (e.g.collectri.h5ad)> --save_dir=<e.g.output/> --dataset=<e.g. replogle> --build_images=<true or false. true for the first time running> --run_test=<true or false. true to run on test data>
```

example command:
```bash
bash scripts/run_grn_evaluation.sh --prediction=resources/grn_models/op/collectri.h5ad --save_dir=output/ --dataset=op --build_images=true --test_run=false
```
