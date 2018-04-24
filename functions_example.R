# load functions
source('functions.R')

###########
# example for using function to test for differences in correlation between two groups 
###########

# data is a matrix of normalized gene expression values (data DO NOT need to be normalized within each group of interest)
# individuals are in columns, genes are in rows
data=read.delim('example_exp_dataset.txt')
# class assignment for each individual, coded as 0 or 1
# individuals should be in the same order in the predictor file and the data matrix
pred=read.delim('example_pred_dataset.txt')

diff_correlation_2grps(data,pred$x,'healthy','sick',plot=T)
