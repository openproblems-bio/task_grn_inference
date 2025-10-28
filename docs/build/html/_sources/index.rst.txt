Documentation for Gene Regulatory Network Inference Benchmark (geneRNIB)
========================================================================


geneRNIB is a living benchmark platform for GRN inference. This platform provides curated datasets for GRN inference and evaluation, standardized evaluation protocols and metrics, computational infrastructure, and a dynamically updated leaderboard to track state-of-the-art methods. 
It runs novel GRNs in the cloud, offers competition scores, and stores them for future comparisons, reflecting new developments over time.

The platform supports the integration of new inference methods, datasets, and protocols. When a new feature is added, previously evaluated GRNs are re-assessed, and the leaderboard is updated accordingly. 
It is designed for both single-modality and multi-omics GRN inference.

.. image:: images/overview.png
   :width: 70%
   :align: center
----

This documentation is supplementary to the paper `geneRNIB: a living benchmark for gene regulatory network inference <https://www.biorxiv.org/content/10.1101/2025.02.25.640181v1.full.pdf>`_ and the `GitHub page <https://github.com/openproblems-bio/task_grn_inference>`_ on the OpenProblems platform. 

- To install geneRNIB, see the `GitHub page <https://github.com/openproblems-bio/task_grn_inference>`_
- To download, see :doc:`dataset` page
- To perform GRN inference using our integrated methods, see :doc:`inference` page
- To run evaluation metrics, see :doc:`evaluation` page
- To extend geneRNIB with new methods, metrics, or datasets, see :doc:`extending` page
- To view the leaderboard of integrated methods, see :doc:`leaderboard` page

.. .. image:: images/grn_models.png
..    :width: 70%
..    :align: center
.. ----


.. Pls see the GitHub page for the list of currently integrated methods. The methods are implemented in Python and R, and they can be used to infer GRNs from the datasets provided by geneRNIB.

.. In addition, three baseline methods are integrated into geneRNIB. These methods are used to evaluate the performance of new methods. The baseline methods are:

.. - **Negative control**: Randomly assigns weights to edges. GRN inference methods should outperform this method.
.. - **Pearson correlation**: Assigns weights based on the Pearson correlation between genes.
.. - **Positive control**: Similar to Pearson correlation with the exception that it uses both inference and evaluation dataset to infer the GRN. This method is expected to outperform most methods.


.. .. list-table:: Authors & contributors
..    :widths: 25 25
..    :header-rows: 1

..    * - name
..      - roles
..    * - Jalil Nourisa
..      - author
..    * - Antoine Passemiers
..      - author
..    * - Robrecht Cannoodt
..      - author
..    * - Marco Stock
..      - author



Contents
--------

.. toctree::
   dataset
   inference
   evaluation
   extending
   leaderboard
  
   