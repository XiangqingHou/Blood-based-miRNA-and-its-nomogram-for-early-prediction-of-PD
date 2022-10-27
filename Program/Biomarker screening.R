library(geepack)
setwd("//Users//houxiangqing//Desktop//temp//Datasets_PPMI//")
#import miRNA
# Reads per million mapped to miRNA normalized read counts for all miRNAs of miRBase v22
rna<-read.csv('mirna_quantification_matrix_rpmmm_norm.txt',header = TRUE,sep='\t')
row.names(rna)<-rna[,1]
rna<-rna[,-1]
rna1<-rna[which(rowSums(rna) > 0),]
#import clinical data
cli<-read.csv('clinical.csv',header = T,sep=',')#training datasets
row.names(cli)<-cli[,5]
cli1<-cli[,-5]
#select target miRNA base on clinical index
rna2<-subset(rna1,select=c(which(colnames(rna1) %in% row.names(cli1))))
rna2<-rna2[which(rowSums(rna2) > 0),]
#rank transfer
rna3<-t(rna2)
#rm expressive level of miRNA equal to 0 more than 50%
t_inData<-t(rna3)
n0 <- apply(t_inData == 0, 1, sum)
i0 <- which(n0 > 1105)
result<-t_inData[-i0,]
result1<-t(result)
#merge miRNA and clinica data
rna4<-merge(result1,cli1,by='row.names',all=TRUE)
row.names(rna4)<-rna4[,1]
rna4<-rna4[,-1]
rna4<-rna4[order(rna4[,'Patient_number']),]
#GEE model to select df gene(HY stage as outcome)
HX_miRNA_Main <- data.frame()
for(k in 1:1081){
  miRNA <- colnames(rna4)[k]
  fit_miRNA <- try(geese(stage ~ rna4[,k] + time + rna4[,k]*time, id = rna4[,'Patient_number'],
                         data=rna4, corstr="exch", family=poisson))
  P_miRNA <- summary(fit_miRNA)[[1]][2,4]
  beta_miRNA <- summary(fit_miRNA)[[1]][2,1]
  HX_miRNA_Main <- rbind(HX_miRNA_Main,data.frame(Gene=miRNA,beta_miRNA=beta_miRNA,P_miRNA=P_miRNA,
                                                  stringsAsFactors=F))
}
#select of dif gene according to P<0.05
diff_miRNA<-HX_miRNA_Main[which(HX_miRNA_Main$P_miRNA<=0.05),]

#GEE model with continuous variable--train
HX_miRNA_Main2<- data.frame()
for(k in 1:1081){
  miRNA <- colnames(rna4)[k]
  fit_miRNA1<-try(geeglm(score_t~rna4[,k] + time + rna4[,k]*time,id=rna4[,'Patient_number'],corstr = "exch",
                         family="gaussian",data=rna4,std.err = "san.se"))
  P_miRNA <- summary(fit_miRNA1)$coefficients[2,4]
  beta_miRNA <- summary(fit_miRNA1)$coefficients[2,1]
  HX_miRNA_Main2 <- rbind(HX_miRNA_Main2,data.frame(Gene=miRNA,beta_miRNA=beta_miRNA,P_miRNA=P_miRNA,
                                                  stringsAsFactors=F))
}
#select of dif gene according to P<0.001
diff_miRNA2<-HX_miRNA_Main2[which(HX_miRNA_Main2$P_miRNA<=0.05),]

#inData intersect model1 and model2
inData<-subset(rna3,select=c(which(colnames(rna3) %in% intersect(diff_miRNA2$Gene,diff_miRNA$Gene))))
#merge diffRNA and clinica data
inData<-inData[order(row.names(inData)),]
cli1<-cli1[order(row.names(cli1)),]
inData1<-merge(inData,cli1,by='row.names',all=TRUE)
row.names(inData1)<-inData1[,1]
inData1<-inData1[,-1]
inData1<-inData1[order(inData1[,'Patient_number']),]
#output datasets
write.csv(inData,file='diff_miRNA.csv',quote=TRUE,row.names=TRUE)