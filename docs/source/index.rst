Documentation for gene regulatory network inference benchmark (geneRNIB) 
==================================================================================================

geneRNIB is a living benchmark platform for GRN inference. This platform provides curated datasets for GRN inference and evaluation, standardized evaluation protocols and metrics, computational infrastructure, and a dynamically updated leaderboard to track state-of-the-art methods. It runs novel GRNs in the cloud, offers competition scores, and stores them for future comparisons, reflecting new developments over time.

The platform supports the integration of new inference methods, datasets and protocols. When a new feature is added, previously evaluated GRNs are re-assessed, and the leaderboard is updated accordingly. The aim is to evaluate both the accuracy and completeness of inferred GRNs. It is designed for both single-modality and multi-omics GRN inference.

In the current version, geneRNIB contains 11 inference methods including both single and multi-omics, 8 evalation metrics, and five datasets (OPSCA, Nakatake, Norman, Adamson, and Replogle).


This documentation is summlementary to the `task_gen_benchmark <https://github.com/openproblems-bio/task_grn_inference>`_ on the OpenProblems platform. You can find out the latest integrated GRN inference methods on this page.

The overview of the benchmark can be found :doc:`overview`.

To integrate your method into geneRNIB, follow the :doc:`extending` section.
Follow the :doc:`inference` section to access inference datasets and infer gene regulatory networks using your method without integrating it into geneRNIB.


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
   inference
   