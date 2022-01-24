# GP_vs_PP

R functions to implement the comparative analysis of genomic and phenomic predictions of growth-related traits in three-way coffee hybrids, as described in [(A. Mbebi et al. 2022)]().

1. The folder code contains the following R scripts:

  * L21_featselect.R which uses the model proposed in [(Nie et al. 2010)](http://papers.nips.cc/paper/3988-efficient-and-robust-feature-selection-via-joint-l21-norms-minimization) to implement GS under the ![equation](https://latex.codecogs.com/gif.latex?%5Ctext%7BL%7D_%7B21%7D)-norm regularized multivariate regression 

  * CV_L21_featselect.R which selects the tuning parameter for L21_featselect.R using K-folds cross-validation

  * Ridge_estim.R compute the Ridge estimate

  * CV_Ridge.R selects the tuning parameter for Ridge_estim.R using K-folds cross-validation

  * L21_joint_estim.R performs GS using the ![equation](https://latex.codecogs.com/gif.latex?%5Ctext%7BL%7D_%7B21%7D)-norm regularized multivariate regression that jointly estimates the regression coefficients and precision matrix
 
  * CV_L21_joint_estim.R selects the tuning parameters for L21_joint_estim.R using K-folds cross-validation

  * MOR.R runs multiple output regression [(He et al. 2016)](https://academic.oup.com/bioinformatics/article/32/12/i37/2288681)

  * CV_MOR.R selects the tuning parameters for MOR.R using K-folds CV

2. The folder data contains data sets for the H1xG and H1xET47 families used in the comparative analysis. 

  * H1xET47_coffeetraits_Step234.csv and H1xG_coffeetraits_Step234.csv contain leaf count (LC), tree height (TH) and trunk diameter (TD) measurements for step2 (i.e. 2 months growth under full and altitude level 600m), step3 (i.e. 2.5 months growth under shade and altitude level 600m) and step4 (i.e. 2 months growth under full sun and altitude level 1300m) for H1xET47 and H1xG families respectively.
  * H1xET47_coffeetraits_Step234.csv
  * H1xET47_coffeetraits_Step234.csv
  * H1xET47_coffeetraits_Step234.csv
  * H1xET47_coffeetraits_Step234.csv

3. The file Example_run.R is an implementation example of all models/functions previously discussed, using simulated data set. 

5. Notes
  * The MTGS.mlasso function is called from [MTGS R package](https://CRAN.R-project.org/package=MTGS) and compute the GEBVs using multivariate LASSO 

  * The MTGS.kmlasso function is called from [MTGS R package](https://CRAN.R-project.org/package=MTGS) and compute the GEBVs using Kernelized multivariate LASSO

  * The MTGS.mrce function is called from [MTGS R package](https://CRAN.R-project.org/package=MTGS) and compute the GEBVs using MRCE

  * By default, all cross-validation scripts use 5-folds

  * Although the codes here were tested on Fedora 29 (Workstation Edition) using R (version 3.6.1), they can run under any Linux or Windows OS distributions, as long as all the required packages are compatible with the desired R version.

  * The following abbreviations are used, CV: Cross Validation, GBLUP: Genomic Best Linear Unbiased Prediction, GS: Genomic Selection, LASSO: Least Absolute Shrinkage and Selection Operator, MOR: Multiple Output Regression
MRCE: Multivariate Regression with Covariance Estimation, MTGS: Multi Traits Genomic Selection.



