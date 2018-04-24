# README

The code provided here implements the method described in [paper info]. 

The goal of the method is to test whether a predictor variable (e.g., environment, disease status, or genotype) affects the degree to which two variables are correlated.

Several key files are provided:

1) example.R is an Rscript that will simulate gene expression data for two groups of individuals, those that are healthy and those that are sick. The pairwise correlation coefficient is decreased in sick individuals relative to healthy individuals for a subset of genes pairs. The script walks through the process of data normalization and implementing a test for differential correlation.

2) functions.R is an Rscript with functions that can be loaded and applied to external datasets.

3) example_functions.R is an Rscript that shows how to use the provided functions.

Please contact Amanda Lea (amandalea7180@gmail.com) with any questions.
