CILP_withLM=function(input_matrix,predictor_DF,predictor_column,covariates_column_start,covariates_column_end,class1_name,class2_name,plot=T) {

# get all possible pairwise combinations of genes
genes=dim(input_matrix)[1]
e_all<-input_matrix
predictor<-predictor_DF[,predictor_column]
pairs<-as.data.frame(t(combn( 1:genes, 2)))

print('scaling data within the classes to be compared')
tmp<-matrix(nrow=dim(e_all)[1],ncol=dim(e_all)[2])
for (i in c(1:dim(tmp)[1])) {
tmp[i,which(predictor==1)]<-scale( t(e_all[i,which(predictor==1)] ))
tmp[i,which(predictor==0)]<-scale( t(e_all[i,which(predictor==0)] ))
} 
e_all_norm<-tmp

print('testing for differences in correlation')

model<-model.matrix(~predictor_DF[,predictor_column] + as.matrix(predictor_DF[,covariates_column_start:covariates_column_end]) )

results <- as.data.frame(t(Reduce(cbind,lapply(1:dim(pairs)[1], function(i) {
tmp1<-(e_all_norm[pairs$V1[i],])
tmp2<-(e_all_norm[pairs$V2[i],])
return( c( cor.test( e_all_norm[pairs$V1[i],which(predictor==0)] , e_all_norm[pairs$V2[i],which(predictor==0)] ,method='spearman')$estimate, cor.test( e_all_norm[pairs$V1[i],which(predictor==1)] , e_all_norm[pairs$V2[i],which(predictor==1)] ,method='spearman')$estimate, summary(lm( scale(tmp1*tmp2) ~ model))$coefficients[2,4]))
}))))
names(results)<-c(paste('cor',class1_name,sep='_'),paste('cor',class2_name,sep='_'),'p_value')
print('saving results')
write.table(results,'results.txt',row.names=F,quote=F,sep='\t')

if(plot==T) {
print('plotting results')
par(mfrow=c(2,2))

# qq-plot
qqplot(-log10(runif(dim(pairs)[1])),-log10(results$p_value),pch=20,xlab='p-value, uniform distribution',ylab='p-value, differential correlation test',bty='n',xlim=c(0,max(-log10(results$p_value))),ylim=c(0,max(-log10(results$p_value))))
x=c(0,1);y=c(0,1);abline(lm(y~x),lty=2)

} 
print('done')
}
