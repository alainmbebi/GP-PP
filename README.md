# GP_vs_PP

R functions to implement the comparative analysis of genomic and phenomic predictions of growth-related traits in three-way coffee hybrids, as described in [( Mbebi A. et al. 2022)]().

1. The folder code contains the R script H1xET47_GP_step234.R in which one can find function for parameters estimation and K-folds cross-validation for all seven statistcal models used in the comparative analysis. This is provided as an example for H1xET47 family only because the results for H1xG family can be obtained by replacing traits and SNP data for H1xET47 with those of H1xG. The same script is used for phenomic prediction, using the Chlorophyll a fluorescence measurements for the corresponding population instead of SNP (except for GBLUP that is not considered for PP).

2. The folder data contains data sets for the H1xG and H1xET47 families used in the comparative analysis. 

  * H1xET47_coffeetraits_Step234.csv and H1xG_coffeetraits_Step234.csv contain measurements for the analysed coffee traits (i.e. LC, TH and TD) over step2, step3 and step4 for H1xET47 and H1xG families respectively.
 
  * H1xET47_Phenomic_Step234.csv and H1xG_Phenomic_Step234.csv contain phenomic (chlorophyll a fluorescence) measurements over step2, step3 and step4 for H1xET47 and H1xG families respectively.
  
  * FinalSNP_H1xET47_impute_mean.csv and FinalSNP_H1xG_impute_mean.csv are SNP data for H1xET47 and H1xG families respectively.
  
  * nugen_spet-gatk_haplotypecaller_SNP.vcf.gz.tbi file that is included for completeness, contains raw SNP data for the two families.

3. Notes

  * By default, all cross-validation scripts use 3-folds 

  * Although the codes here were tested on Fedora 29 (Workstation Edition) using R (version 3.6.1), they can run under any Linux or Windows OS distributions, as long as all the required packages are compatible with the desired R version.

  * The following abbreviations are used, LC: leaf count, TH: tree height, TD: trunk diameter, GP: genomic prediction. PP: phenomic prediction, L21-joint:  equation-norm regularized multivariate regression, RR: ridge regression , mLASSO: multiple LASSO, EN: elastic-net, BL: Bayesian LASSO, GBLUP: Genomic Best Linear Unbiased Prediction and mBayesB: multiple-trait BayesB.
  
  * Step2, step3 and step4 correspond respectively to: 2 months growth under full sun altitude level 600m, 2.5 months growth under shade altitude level 600m and 2 months growth under full sun and altitude level 1300m.



