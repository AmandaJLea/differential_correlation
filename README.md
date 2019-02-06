# README

The code provided here implements the method described in "Genetic and environmental perturbations lead to regulatory decoherence" (https://www.biorxiv.org/content/early/2018/07/14/369306). The method is referred to as "Correlation by Individual Level Product" or CILP.

The goal of the method is to test whether a predictor variable (e.g., environment, disease status, or genotype) affects the degree to which two variables are correlated. 

Our approach is based on the fact that the Pearson correlation coefficient is equal to the average element-wise product of two traits measured across individuals, after each trait is mean centered and scaled; this value reflects the relationship between two traits within a population sample. By extension, to obtain a measure of the degree of correlation between two traits for each individual in a sample, we  simply keep the vector of products and do not perform averaging across individuals. This vector of products can be used as the outcome variable for a variety of statistical tests (e.g., linear model or linear mixed model). 

Note that in the manuscript referenced above, we describe two versions of CILP. The main version is used to test for effects of infection on pairwise gene expression correlations and for effects of metabolic syndrome on pairwise metabolite correlations, while a modified version is used for scalibility during correlation QTL screening. The code provided here focuses on the main version of CILP.  

Several files are provided:

1) example.R is an Rscript that will simulate gene expression data for two groups of individuals, those that are healthy and those that are sick. The pairwise correlation coefficient is decreased in sick individuals relative to healthy individuals for a subset of genes pairs. The script walks through the process of data normalization, implementing CILP, and checking null expectations. It is strongly reccomended that any new user walks through this script first, before attempting to implement CILP on their own data.

2) functions.R is an Rscript with functions that can be loaded and applied to external datasets.

3) functions_example.R is an Rscript that shows how to use the provided functions.

Please contact Amanda Lea (amandalea7180@gmail.com) with any questions.
