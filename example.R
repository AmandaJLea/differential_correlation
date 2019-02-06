############
# set parameters and underlying correlation structure
# in this example, correlations are reduced in magnitude in sick individuals relative to healthy individuals for a subset of gene pairs
############

n=1000
healthy=n/2
sick=n/2
prop_change=1/10
change_effect_size=1/4
genes=100

mat_h<-diag(genes)
mat_s<-diag(genes)

# gene pairs in healthy individuals display a range of correlation values
tmp<- sample( runif(10000,0.1,0.9) ,length(mat_h[upper.tri(mat_h)]))
mat_h[upper.tri(mat_h)] <- tmp
mat_h[lower.tri(mat_h)] <- t(mat_h)[lower.tri(mat_h)]

# for sick people, correlations are similar but reduced at a subset of sites
tmp<-as.data.frame(tmp)
tmp$change<-0
tmp$change[sample(1:dim(tmp)[1],dim(tmp)[1]*prop_change)]<-1
tmp$new<-tmp$tmp
tmp$new[which(tmp$change==1)]<-tmp$tmp[which(tmp$change==1)] * change_effect_size

mat_s[upper.tri(mat_s)] <- tmp$new
mat_s[lower.tri(mat_s)] <- t(mat_s)[lower.tri(mat_s)]

############
# simulate gene expression data from underlying correlation structure
############

library(MASS)
library(Matrix)
mat_h2<-as.matrix(nearPD(mat_h,ensureSymmetry=TRUE)[[1]])
mat_s2<-as.matrix(nearPD(mat_s,ensureSymmetry=TRUE)[[1]])

# simulate mean differences in expression between healthy and sick individuals for a subset of genes as well (to make things more realistic, and to check normalization)
mean_h<-rep(0,genes)
mean_s<-rep(0,genes)
mean_s[sample(1:genes,genes*prop_change)]<-0.25

e_h = mvrnorm(healthy, mean_h, mat_h2) 
e_s = mvrnorm(sick, mean_s, mat_s2) 
e_all<-as.data.frame(cbind(t(e_h),t(e_s)))

# confirm mean differences in expression
predictor<-c(rep(0,healthy),rep(1,sick))
p1<-apply(e_all,1,function(x) summary(lm(x ~ predictor))$coefficients[2,4])
qqplot(-log10(runif(genes)),-log10(p1),pch=20,xlab='p-value, uniform distribution',ylab='p-value, mean effect',bty='n',xlim=c(0,max(-log10(p1))),ylim=c(0,max(-log10(p1))))
x=c(0,1);y=c(0,1);abline(lm(y~x),lty=2)

############
# test for differences in correlation between healthy and sick
############

# get all possible pairwise combinations of genes
pairs<-as.data.frame(t(combn( 1:genes, 2)))

# scale data within the classes to be compared
tmp<-matrix(nrow=dim(e_all)[1],ncol=dim(e_all)[2])
for (i in c(1:dim(tmp)[1])) {
tmp[i,which(predictor==1)]<-scale( t(e_all[i,which(predictor==1)] ),center = TRUE, scale = TRUE)
tmp[i,which(predictor==0)]<-scale( t(e_all[i,which(predictor==0)] ),center = TRUE, scale = TRUE)
} 
e_all_norm<-tmp

# confirm there are no longer mean differences in expression
p2<-apply(e_all_norm,1,function(x) summary(lm(x ~ predictor))$coefficients[2,4])
summary(p2)

# confirm there are no/minimal variance differences in expression either
p2<-apply(e_all_norm,1,function(x) var(x[which(predictor==0)]) - var(x[which(predictor==1)]))
summary(p2)

# test for differences in correlation (this will take a few minutes to run, depending on sample size and gene number)
results <- as.data.frame(t(Reduce(cbind,lapply(1:dim(pairs)[1], function(i) {
tmp1<-(e_all_norm[pairs$V1[i],])
tmp2<-(e_all_norm[pairs$V2[i],])
# here, a linear model is used to test for differences in correlation between groups; however, any number of other statistical methods (e.g., linear models with additional covariates or linear mixed effects models) could be used in the same framework
# to use an alternate statistical approach, extract the product matrix and use it as the outcome variable in another type of statistical model  
return( c( cor.test( e_all_norm[pairs$V1[i],which(predictor==0)] , e_all_norm[pairs$V2[i],which(predictor==0)] ,method='spearman')$estimate, cor.test( e_all_norm[pairs$V1[i],which(predictor==1)] , e_all_norm[pairs$V2[i],which(predictor==1)] ,method='spearman')$estimate, summary(lm( scale(tmp1*tmp2) ~ predictor))$coefficients[2,4]))
}))))
names(results)<-c('cor_healthy','cor_sick','p_value')

############
# visualize results
############

par(mfrow=c(2,2))

# volcano plot
# remember we simulated a loss of correlation in sick individuals, so this plot should be lopsided
library(qvalue)
results$qvalue<-qvalue(results$p_value)$qvalues
plot(results$cor_healthy - results$cor_sick,-log10(results$p_value),bty='n',col=as.factor(results$qvalue<0.05),xlab='difference in Spearman correlation (healthy-sick)',ylab='-log10, p-value')

# qq-plot
qqplot(-log10(runif(dim(pairs)[1])),-log10(results$p_value),pch=20,xlab='p-value, uniform distribution',ylab='p-value, differential correlation test',bty='n',xlim=c(0,max(-log10(results$p_value))),ylim=c(0,max(-log10(results$p_value))))
x=c(0,1);y=c(0,1);abline(lm(y~x),lty=2)

# example
i<-which( results$p_value==min(results$p_value))
plot( e_all_norm[pairs$V1[i],which(predictor==0)] , e_all_norm[pairs$V2[i],which(predictor==0)] ,bty='n',xlab='Normalized expression of gene 1',ylab='Normalized expression of gene 2',main='Healthy group')
abline(lm ( e_all_norm[pairs$V2[i],which(predictor==0)] ~ e_all_norm[pairs$V1[i],which(predictor==0)] ),lty=2)

plot( e_all_norm[pairs$V1[i],which(predictor==1)] , e_all_norm[pairs$V2[i],which(predictor==1)] ,bty='n',xlab='Normalized expression of gene 1',ylab='Normalized expression of gene 2',main='Sick group')
abline(lm ( e_all_norm[pairs$V2[i],which(predictor==1)] ~ e_all_norm[pairs$V1[i],which(predictor==1)] ),lty=2)

############
# check null expectations through permutation
# permute the sample labels before manipulating the gene expression data
############

predictor_perm<-sample(predictor)

# scale data within the classes to be (permuted) compared
tmp<-matrix(nrow=dim(e_all)[1],ncol=dim(e_all)[2])
for (i in c(1:dim(tmp)[1])) {
tmp[i,which(predictor_perm==1)]<-scale( t(e_all[i,which(predictor_perm==1)] ),center = TRUE, scale = TRUE)
tmp[i,which(predictor_perm==0)]<-scale( t(e_all[i,which(predictor_perm==0)] ),center = TRUE, scale = TRUE)
} 
e_all_norm<-tmp

# test for differences in correlation using permuted predictor variable (this will take a few minutes to run, depending on sample size and gene number)
results_perm <- as.data.frame(t(Reduce(cbind,lapply(1:dim(pairs)[1], function(i) {
tmp1<-(e_all_norm[pairs$V1[i],])
tmp2<-(e_all_norm[pairs$V2[i],])
return( c( cor.test( e_all_norm[pairs$V1[i],which(predictor_perm==0)] , e_all_norm[pairs$V2[i],which(predictor_perm==0)] ,method='spearman')$estimate, cor.test( e_all_norm[pairs$V1[i],which(predictor_perm==1)] , e_all_norm[pairs$V2[i],which(predictor_perm==1)] ,method='spearman')$estimate, summary(lm( scale(tmp1*tmp2) ~ predictor_perm))$coefficients[2,4]))
}))))
names(results_perm)<-c('cor_healthy','cor_sick','p_value')

# qq-plot - should not show any enrichment of low p-values relative to a uniform distribution  
qqplot(-log10(runif(dim(pairs)[1])),-log10(results_perm$p_value),pch=20,xlab='p-value, uniform distribution',ylab='p-value, differential correlation test (permuted)',bty='n',xlim=c(0,max(-log10(results_perm$p_value))),ylim=c(0,max(-log10(results_perm$p_value))))
x=c(0,1);y=c(0,1);abline(lm(y~x),lty=2)

############
# check null expectations through permutation
# permute the sample labels before running statistical tests (shown here = linear models)
############

predictor_perm<-sample(predictor)

# scale data within the (non permuted) classes to be compared
tmp<-matrix(nrow=dim(e_all)[1],ncol=dim(e_all)[2])
for (i in c(1:dim(tmp)[1])) {
tmp[i,which(predictor==1)]<-scale( t(e_all[i,which(predictor==1)] ),center = TRUE, scale = TRUE)
tmp[i,which(predictor==0)]<-scale( t(e_all[i,which(predictor==0)] ),center = TRUE, scale = TRUE)
} 
e_all_norm<-tmp

# test for differences in correlation using permuted predictor variable (this will take a few minutes to run, depending on sample size and gene number)
results_perm2 <- as.data.frame(t(Reduce(cbind,lapply(1:dim(pairs)[1], function(i) {
tmp1<-(e_all_norm[pairs$V1[i],])
tmp2<-(e_all_norm[pairs$V2[i],])
return( c( cor.test( e_all_norm[pairs$V1[i],which(predictor_perm==0)] , e_all_norm[pairs$V2[i],which(predictor_perm==0)] ,method='spearman')$estimate, cor.test( e_all_norm[pairs$V1[i],which(predictor_perm==1)] , e_all_norm[pairs$V2[i],which(predictor_perm==1)] ,method='spearman')$estimate, summary(lm( scale(tmp1*tmp2) ~ predictor_perm))$coefficients[2,4]))
}))))
names(results_perm2)<-c('cor_healthy','cor_sick','p_value')

# qq-plot - should not show any enrichment of low p-values relative to a uniform distribution  
qqplot(-log10(runif(dim(pairs)[1])),-log10(results_perm2$p_value),pch=20,xlab='p-value, uniform distribution',ylab='p-value, differential correlation test (permuted)',bty='n',xlim=c(0,max(-log10(results_perm2$p_value))),ylim=c(0,max(-log10(results_perm2$p_value))))
x=c(0,1);y=c(0,1);abline(lm(y~x),lty=2)
