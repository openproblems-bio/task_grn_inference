Welcome to the documentation for the GRN inference task on OpenProblems Benchmarking platform!
===================================

This task is one of many other tasks hosted on the `OpenProblems benchmarking platform <https://openproblems.bio/>`.
The interactive benchmarking results are hosted `here <https://openproblems.bio/results/>`.

The *GRN inference* task focuses on the inference of gene regulatory networks (GRN) from RNA-Seq expression or chromatin accessibility data (ATAC-Seq) or both. 
The pipeline evaluates inferred GRNs against pertubation data, by training two types of regression models. This type of evaluation is closer to evaluating the biological knowledge that a GRN should represent, instead of evaluating the presence of edges in a statistical way only, as commonly done by using metrics, such as AUPRC or AUROC.

Jump to the :doc:`evaluation` section to get a deeper explanation on how the benchmarking task is setup.
If you want to add your own datasets or algorithms to the benchmark, check our the :doc:`add_things` section.

.. note::

   This project is under active development and this documentation is still a draft.

Contents
--------

.. toctree::

   evaluation
   add_stuff