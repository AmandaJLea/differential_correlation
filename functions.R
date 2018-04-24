diff_correlation_2grps=function(input_matrix,input_predictor,class1_name,class2_name,plot=T) {

# get all possible pairwise combinations of genes
genes=dim(data)[1]
e_all<-input_matrix
predictor<-input_predictor
pairs<-as.data.frame(t(combn( 1:genes, 2)))

print('scaling data within the classes to be compared')
tmp<-matrix(nrow=dim(e_all)[1],ncol=dim(e_all)[2])
for (i in c(1:dim(tmp)[1])) {
tmp[i,which(predictor==1)]<-scale( t(e_all[i,which(predictor==1)] ))
tmp[i,which(predictor==0)]<-scale( t(e_all[i,which(predictor==0)] ))
} 
e_all_norm<-tmp

print('testing for differences in correlation')
results <- as.data.frame(t(Reduce(cbind,lapply(1:dim(pairs)[1], function(i) {
tmp1<-(e_all_norm[pairs$V1[i],])
tmp2<-(e_all_norm[pairs$V2[i],])
# here, a wilcoxon test is used to test for differences in correlation between groups; however, any number of other statistical methods (e.g., linear models or linear mixed effects models) could be used in the same framework
return( c( cor.test( e_all_norm[pairs$V1[i],which(predictor==0)] , e_all_norm[pairs$V2[i],which(predictor==0)] ,method='spearman')$estimate, cor.test( e_all_norm[pairs$V1[i],which(predictor==1)] , e_all_norm[pairs$V2[i],which(predictor==1)] ,method='spearman')$estimate, wilcox.test( (tmp1*tmp2) ~ predictor)$p.value))
}))))
names(results)<-c('cor_healthy','cor_sick','p_value')
print('saving results')
write.table(results,'results.txt',row.names=F,quote=F,sep='\t')

if(plot==T) {
print('plotting results')
par(mfrow=c(2,2))

# volcano plot
library(qvalue)
results$qvalue<-qvalue(results$p_value)$qvalues
plot(results$cor_healthy - results$cor_sick,-log10(results$p_value),bty='n',col=as.factor(results$qvalue<0.05),xlab='difference in correlation',ylab='-log10, p-value')

# qq-plot
qqplot(-log10(runif(dim(pairs)[1])),-log10(results$p_value),pch=20,xlab='p-value, uniform distribution',ylab='p-value, differential correlation test',bty='n',xlim=c(0,max(-log10(results$p_value))),ylim=c(0,max(-log10(results$p_value))))
x=c(0,1);y=c(0,1);abline(lm(y~x),lty=2)

# example
i<-which( results$p_value==min(results$p_value))
plot( e_all_norm[pairs$V1[i],which(predictor==0)] , e_all_norm[pairs$V2[i],which(predictor==0)] ,bty='n',xlab='Normalized expression of gene 1',ylab='Normalized expression of gene 2',main=class1_name)
abline(lm ( e_all_norm[pairs$V2[i],which(predictor==0)] ~ e_all_norm[pairs$V1[i],which(predictor==0)] ),lty=2)

plot( e_all_norm[pairs$V1[i],which(predictor==1)] , e_all_norm[pairs$V2[i],which(predictor==1)] ,bty='n',xlab='Normalized expression of gene 1',ylab='Normalized expression of gene 2',main=class2_name)
abline(lm ( e_all_norm[pairs$V2[i],which(predictor==1)] ~ e_all_norm[pairs$V1[i],which(predictor==1)] ),lty=2)

} 
print('done')
}
  
