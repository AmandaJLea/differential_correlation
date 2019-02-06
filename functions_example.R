# load functions
source('functions.R')

###########
# example for using function to test for differences in correlation between two groups 
# this version expects covariates
# a linear model is used to test for differences in correlation as a function of the predictor variable, while controlling for covariates
###########

# data is a matrix of normalized gene expression values (data DO NOT need to be normalized within each group of interest; however, expression values for each gene should approximate a normal distribution, they should not be in raw count format)
# individuals are in columns, genes are in rows
data=read.delim('example_exp_dataset.txt')

# class assignment for each individual, coded as 0 or 1, is in the first column
# individuals should be in the same order in the predictor file and the data matrix
yvars=read.delim('example_pred_dataset.txt')

# covariates are also expected in additional columns in the predictor file

# CILP_withLM=function(input_matrix,predictor_DF,predictor_column,covariates_column_start,covariates_column_end,class1_name,class2_name,plot=T) 

CILP_withLM(data,yvars,1,2,3,'healthy','sick',plot=T)

###########
# example for using function to test for differences in correlation between two groups 
# this version DOES NOT expect covariates
# a linear model is used to test for differences in correlation as a function of the predictor variable
###########

data=read.delim('example_exp_dataset.txt')
yvars=read.delim('example_pred_dataset.txt')

# CILP_withLM_noCOV=function(input_matrix,predictor_DF,predictor_column,class1_name,class2_name,plot=T) 

CILP_withLM_noCOV(data,yvars,1,'healthy','sick',plot=T)

###########
# example for using function to plot differences in correlation for a pair of genes of interest
###########

data=read.delim('example_exp_dataset.txt')
yvars=read.delim('example_pred_dataset.txt')

# output from running functions described above
results=read.delim('results.txt')
pairs=read.delim('pairs_tested.txt')

# chose an interesting candidate gene pair
library(qvalue)
pairs[which(qvalue(results$p_value)$qvalues<0.05)[1],]

# plot_pair=function(input_matrix,predictor_DF,predictor_column,gene1,gene2,class1_name,class2_name)
plot_pair(data,yvars,1,6,72,'healthy','sick')
