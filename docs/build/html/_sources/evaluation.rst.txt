
GRN evaluation
=================
The evaluation metrics used in geneRNIB are summarized below. 

  
.. image:: images/metrics.png
   :width: 90%
   :align: center
----

Each metric has particular implementation requirements that limit its applicability to certain datasets. The datasets applicable to each metric are listed in the following table.

In addition, to ensure that only informative metrics contribute to the benchmark, each metric is evaluated per dataset against two quality control criteria:

1. **Performance threshold** — the best score across all methods must exceed a threshold, defined as the maximum of a global floor (metric-specific) and the negative control score for that dataset. This filters out metrics where no method achieves meaningful signal.
2. **Variability** — the coefficient of variation (CV = range / mean) across methods must be ≥ 0.2. This filters out metrics that fail to discriminate between methods.

A metric is retained for a given dataset only if it passes both criteria. Those that passed quality control are marked with a star (★) in the table below. Note that the quality control process is repeated after adding a new GRN inference method to reflect the updated performance comparison.

.. image:: images/metric_quality_evaluation.png
   :width: 100%
   :align: center
----


For a detailed description of each metric, refer to the geneRNIB paper. To map the naming conventions used in the code and this table, refer to `surrogate_names` in `config.py` file in the `src/utils/` directory.

The evaluation metrics expect the inferred network to be in the form of an AnnData object with specific format as explained here. 
It should be noted that the metric currently evaluate only the **top TF-gene pairs**, currently limited to **50,000 edges**, ranked by their assigned weight.  

The inferred network should have a tabular format with the following columns:  

  - `source`: TF gene name
  - `target`: Target gene gene  
  - `weight`: Regulatory importance/likelihood score/etc.  

See `resources/grn_benchmark/prior/collectri.h5ad` for an example of the expected format.



Running GRN evaluation without docker
----------------------------------------
Considering that Docker is not supported by certtain systems, you can run the evaluation without Docker by following these steps:

```bash
bash src/metrics/all_metrics/run_local.sh --dataset <dataset_name> --prediction=<inferred GRN (e.g.collectri.h5ad)> --score <output_score_file.h5ad> --num_workers <number_of_workers>
```

example command:

```bash
bash src/metrics/all_metrics/run_local.sh --dataset op --prediction=resources/grn_models/op/collectri.h5ad --score=output_score_file.h5ad --num_workers=20
```

If you are evaluating a new GRN model, which is not part of geneRNIB, make you to generate the consensus prior file for the dataset you are evaluating on. 

```bash
bash scripts/prior/run_consensus.sh --dataset op --new_model {new_grn_model_file.h5ad}
```

This will add your model to the previous ones and create the new consensus prior file needed for evaluation.

Running GRN evaluation using standard pipeline
----------------------------------------

To run the evalution for a given GRN and dataset, use the following command:

```bash
bash scripts/run_grn_evaluation.sh --prediction=<inferred GRN (e.g.collectri.h5ad)> --save_dir=<e.g.output/> --dataset=<e.g. replogle> --build_images=<true or false. true for the first time running> 
```

example command:

```bash
bash scripts/run_grn_evaluation.sh --prediction=resources/grn_models/op/collectri.h5ad --save_dir=output/ --dataset=op --build_images=true 
```

