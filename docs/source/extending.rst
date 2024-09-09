Extending the pipeline
======================

Currently the perturbation dataset from the Open Problems: Single Cell Perturbation 2023 data science competition is used.
It provides single-cell perturbation gene expression for peripheral blood mononuclear cells (PBMCs), along with integrated multi-omics data of scRNA-seq and scATAC-seq for the baseline compound from the same experiment.
It includes 146 perturbations, making it the largest drug perturbation study on primary human tissue with donor replicates. 

Currently, the following six enhancer aware GRN inference methods (eGRN methods) are implemented in the pipeline:

#. Scenic+ (`Paper <https://doi.org/10.1038/s41592-023-01938-4>`_)
#. CellOracle (`Paper <https://doi.org/10.1038/s41586-022-05688-9>`_)
#. FigR (`Paper <https://doi.org/10.1016/j.xgen.2022.100166>`_)
#. scGLUE (`Paper <https://doi.org/10.1038/s41587-022-01284-4>`_)
#. GRaNIE (`Paper <https://doi.org/10.15252/msb.202311627>`_)
#. ANANSE (`Paper <https://doi.org/10.1093/nar/gkab598>`_)

To add a method to the repository, follow the instructions in the ``scripts/add_a_method.sh`` script.




