Objects
=======

This provides an overview of the objects used in the GRN inference task.

.. toctree::

   objects

Perturbation
------------

Perturbation dataset for benchmarking.

Example file: ``resources/grn-benchmark/perturbation_data.h5ad``

Format:

.. list-table:: AnnData object
   :header-rows: 1

   * - Attribute
     - Description
   * - ``obs``
     - 'cell_type', 'sm_name', 'donor_id', 'plate_name', 'row', 'well', 'cell_count'
   * - ``layers``
     - 'n_counts', 'pearson', 'lognorm'

Slot description:

.. list-table:: 
   :header-rows: 1

   * - Slot
     - Type
     - Description
   * - ``obs["cell_type"]``
     - ``string``
     - The annotated cell type of each cell based on RNA expression.
   * - ``obs["sm_name"]``
     - ``string``
     - The primary name for the (parent) compound (in a standardized representation) as chosen by LINCS. This is provided to map the data in this experiment to the LINCS Connectivity Map data.
   * - ``obs["donor_id"]``
     - ``string``
     - Donor id.
   * - ``obs["plate_name"]``
     - ``string``
     - Plate name 6 levels.
   * - ``obs["row"]``
     - ``string``
     - Row name on the plate.
   * - ``obs["well"]``
     - ``string``
     - Well name on the plate.
   * - ``obs["cell_count"]``
     - ``string``
     - Number of single cells pseudobulked.
   * - ``layers["n_counts"]``
     - ``double``
     - Pseudobulked values using mean approach.
   * - ``layers["pearson"]``
     - ``double``
     - Normalized values using pearson residuals.
   * - ``layers["lognorm"]``
     - ``double``
     - Normalized values using shifted logarithm.

Control Method
--------------

Path:
``src/methods``

A control method.

Arguments:

.. list-table:: 
   :header-rows: 1

   * - Name
     - Type
     - Description
   * - ``--perturbation_data``
     - ``file``
     - Perturbation dataset for benchmarking.
   * - ``--layer``
     - ``string``
     - Which layer of perturbation data to use to find tf-gene relationships. Default: ``lognorm``.
   * - ``--prior_data``
     - ``file``
     - Prior data used for GRN benchmark.
   * - ``--prediction``
     - ``file``
     - (*Output*) GRN prediction.

Label
-----

Path:
``src/metrics``

A metric to evaluate the performance of the inferred GRN

Arguments:

.. list-table:: 
   :header-rows: 1

   * - Name
     - Type
     - Description
   * - ``--perturbation_data``
     - ``file``
     - Perturbation dataset for benchmarking.
   * - ``--prediction``
     - ``file``
     - GRN prediction.
   * - ``--score``
     - ``file``
     - (*Output*) File indicating the score of a metric.
   * - ``--reg_type``
     - ``string``
     - (*Optional*) name of regression to use. Default: ``ridge``.
   * - ``--subsample``
     - ``integer``
     - (*Optional*) number of samples randomly drawn from perturbation data. Default: ``-1``.


GRN
---

GRN prediction

Example file: ``resources/grn-benchmark/collectri.csv``

Format:

.. list-table:: Tabular data
   :header-rows: 1

   * - Attribute
     - Description
   * - ``Columns``
     - 'source', 'target', 'weight'

Slot description:

.. list-table:: 
   :header-rows: 1

   * - Column
     - Type
     - Description
   * - ``source``
     - ``string``
     - Source of regulation.
   * - ``target``
     - ``string``
     - Target of regulation.
   * - ``weight``
     - ``float``
     - Weight of regulation.


Score
-----

File indicating the score of a metric.

Example file: ``resources/grn-benchmark/score.csv``

Format:

.. list-table:: Tabular data
   :header-rows: 1

   * - Attribute
     - Description
   * - ``Columns``
     - 'accuracy', 'completeness'

Slot description:

.. list-table:: 
   :header-rows: 1

   * - Column
     - Type
     - Description
   * - ``accuracy``
     - ``string``
     - (*Optional*) some explanation.
   * - ``completeness``
     - ``double``
     - (*Optional*) some explanation.


Prior data
----------

Prior data used for grn benchmark

Example file: ``resources/grn-benchmark/prior_data.h5ad``

Format:

.. list-table:: AnnData object
   :header-rows: 1

   * - Attribute
     - Description
   * - ``uns``
     - 'tf_list'

Slot description:

.. list-table:: 
   :header-rows: 1

   * - Slot
     - Type
     - Description
   * - ``uns["tf_list"]``
     - ``list``
     - List of known TFs obtained from https://resources.aertslab.org/cistarget.


Multiomics RNA
--------------

RNA expression for multiomics data.

Example file: ``resources/grn-benchmark/multiomics_rna.h5ad``

Format:

.. list-table:: AnnData object
   :header-rows: 1

   * - Attribute
     - Description
   * - ``obs``
     - 'cell_type', 'donor_id'

Slot description:

.. list-table:: 
   :header-rows: 1

   * - Slot
     - Type
     - Description
   * - ``obs["cell_type"]``
     - ``string``
     - The annotated cell type of each cell based on RNA expression.
   * - ``obs["donor_id"]``
     - ``string``
     - Donor id.

Method
------

Path:
``src/methods``

A GRN inference method

Arguments:

.. list-table:: 
   :header-rows: 1

   * - Name
     - Type
     - Description
   * - ``--multiomics_rna``
     - ``file``
     - RNA expression for multiomics data.
   * - ``--multiomics_atac``
     - ``file``
     - (*Optional*) Peak data for multiomics data.
   * - ``--prediction``
     - ``file``
     - (*Output*) GRN prediction.


Multiomics ATAC
---------------

Peak data for multiomics data.

Example file: ``resources/grn-benchmark/multiomics_atac.h5ad``

Format:

.. list-table:: AnnData object
   :header-rows: 1

   * - Attribute
     - Description
   * - ``obs``
     - 'cell_type', 'donor_id'

Slot description:

.. list-table:: 
   :header-rows: 1

   * - Slot
     - Type
     - Description
   * - ``obs["cell_type"]``
     - ``string``
     - The annotated cell type of each cell based on RNA expression.
   * - ``obs["donor_id"]``
     - ``string``
     - Donor id.









