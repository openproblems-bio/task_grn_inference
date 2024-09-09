Evaluation
==========

The evaluation is done with the help of pertubation data, using two different approaches:

#. Regression from GRN regulations to target expression
#. Regression from TF expression of predicted regulators to target expression

|

.. image:: images/regressions.png
   :width: 100 %
   :alt: overview of the two regression evaluation approaches
   :align: center

|
|


Evaluation 1: Regression from GRN regulations to target expression
------------------------------------------------------------------
The first approach we used is similar to GRaNPA and the multivariate decision tree in Decoupler, where regulatory weights from the GRN form the feature space to predict perturbation data. In this method, we train one model per sample. The feature space matrix has dimensions of genes by transcription factors (TFs), with values being the regulatory weights from the GRN or 0 if the link is absent. The target space matrix represents the perturbation data for each sample. We evaluate the model's predictive performance using a 5-fold cross-validation scheme and the coefficient of determination (R²) as the metric. LightGBM is used for computational efficiency.


Evaluation 2: Regression from TF expression of predicted regulators to target expression
----------------------------------------------------------------------------------------
In the second approach, instead of using regulatory weights, we utilized the expression of putative regulators (TFs) from the perturbation data to construct the feature space. We fit one model per gene, selecting regulators based on the regulatory weights suggested by the GRNs. This method is similar to many modern GRN inference techniques.





