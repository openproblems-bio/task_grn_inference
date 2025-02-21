Documentation for gene regulatory network inference benchmark (geneRNIB) 
==================================================================================================

geneRNIB is a living benchmark platform for GRN inference. This platform provides curated datasets for GRN inference and evaluation, standardized evaluation protocols and metrics, computational infrastructure, and a dynamically updated leaderboard to track state-of-the-art methods. It runs novel GRNs in the cloud, offers competition scores, and stores them for future comparisons, reflecting new developments over time.

The platform supports the integration of new inference methods, datasets and protocols. When a new feature is added, previously evaluated GRNs are re-assessed, and the leaderboard is updated accordingly. The aim is to evaluate both the accuracy and completeness of inferred GRNs. It is designed for both single-modality and multi-omics GRN inference.

In the current version, geneRNIB contains 11 inference methods including both single and multi-omics, 8 evalation metrics, and five datasets (OPSCA, Nakatake, Norman, Adamson, and Replogle).


This documentation is summlementary to the `task_gen_benchmark <https://github.com/openproblems-bio/task_grn_inference>`_ on the OpenProblems platform. You can find out the latest integrated GRN inference methods on this page.
The results of the benchmark can be found `here <https://openproblems.bio/results/>`_.
The paper outlining this benchmark can be found `here <?>`_.

The overview of the benchmark can be found :doc:`overview`.




This task is one of many other tasks hosted on the `OpenProblems benchmarking platform <https://openproblems.bio/>`_.
The interactive benchmarking results are hosted `here <https://openproblems.bio/results/>`_.

The **GRN inference** task focuses on the inference of gene regulatory networks (GRN) from RNA-Seq expression or chromatin accessibility data (ATAC-Seq) or both. 
The pipeline evaluates inferred GRNs against pertubation data, by training two types of regression models. This type of evaluation is closer to evaluating the biological knowledge that a GRN should represent, instead of evaluating the presence of edges in a statistical way only, as commonly done by using metrics, such as AUPRC or AUROC.

Jump to the :doc:`overview` section to get a first summary of the pipeline.
If you want to add your own datasets or algorithms to the benchmark, check our the :doc:`extending` section.


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


.. note::

   This project is under active development and this documentation is still a draft.

Contents
--------

.. toctree::

   overview
   evaluation
   extending
   