Overview
========

The pipeline evaluates inferred GRNs against perturbation data. The evaluation is done by training two types of regression models. 
This type of evaluation is closer to evaluating the biological knowledge that a GRN should represent, instead of evaluating the presence of edges in a statistical way only, as commonly done by using metrics, such as AUPRC or AUROC. 
The general setup is shown in the figure below:

.. image:: images/overview.png
   :width: 100 %
   :alt: overview of pipeline
   :align: center

The pipeline can evaluate algorithms that leverage only one of the multi-omic data types (RNA-Seq or ATAC-Seq) or both.
It also evaluates the performance of two controls:

#. As a *negative control*, the pipeline evaluates the performance of a random network.
#. As a *positive control*, the pipeline evaluates the performance of a network derived from correlation of genes in the perturbation dataset used for evaluation.

The two types of regression models are:

#. Regression from GRN regulations to target expression
#. Regression from TF expression of predicted regulators to target expression

More details can be found in section :doc:`evaluation`.

In the future, other prior knowledge can be incorporated in the pipeline to allow for evaluation with binary classification metrics and as an additional control method. 

Installation
------------

You need to have Docker, Java, and Viash installed. Follow `these instructions <https://openproblems.bio/documentation/fundamentals/requirements/>`_ to install the required dependencies.

Download resources
------------------

.. code-block:: bash

    git clone git@github.com:openproblems-bio/task_grn_benchmark.git

    cd task_grn_benchmark

    # download resources
    scripts/download_resources.sh

Infer a GRN
-----------

.. code-block:: bash

    viash run src/methods/dummy/config.vsh.yaml -- --multiomics_rna resources/grn-benchmark/multiomics_rna.h5ad --multiomics_atac resources/grn-benchmark/multiomics_atac.h5ad --prediction output/dummy.csv


Similarly, run the command for other methods.

Evaluate a GRN
--------------

.. code-block:: bash
    
    scripts/run_evaluation.sh --grn resources/grn-benchmark/grn_models/collectri.csv 

Similarly, run the command for other GRN models.