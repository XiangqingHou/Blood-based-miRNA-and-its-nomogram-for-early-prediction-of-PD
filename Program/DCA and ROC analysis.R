# install dcurves to perform DCA from CRAN
install.packages("dcurves")
# install other packages used in this tutorial
install.packages(
  c("tidyverse", "survival", "gt", "broom",
    "gtsummary", "rsample", "labelled")
)
# load packages
library(dcurves)
library(tidyverse)
library(gtsummary)
library(rms)
################################import data
setwd("//Users//houxiangqing//Desktop//temp//Datasets_PPMI//")
spline <- read.csv("train.csv")
test <- read.csv("External1.csv")
test$Education_Level<-factor(test$Education_Level,labels = c('<12y','12-16y','>16y','Unknown'))
test$Gender<-factor(test$Gender,labels = c("Man","Woman"))
test1<-datadist(test)
options(datadist='test1')
spline$Education_Level<-factor(spline$Education_Level,labels = c('<12y','12-16y','>16y','Unknown'))
spline$Gender<-factor(spline$Gender,labels = c("Man","Woman"))
spline1<-datadist(spline)
options(datadist='spline1')

###################################formulate model
#nomogram model
library(InformationValue)
library(ROCR)
library(arulesViz)
library(pROC)
mod1 <- glm(outcome ~ Gender+Age+Education_Level+Transcriptional_Score, data = spline, family = binomial)
#miRNA model
mod2 <- glm(outcome ~ Transcriptional_Score, data = spline, family = binomial)
#Clinical model
mod3 <- glm(outcome ~ Gender+Age+Education_Level, data = spline, family = binomial)
############################AUROC-test datasets
#nomogram
test.probs=predict(mod1,type="response",newdata=test)
pred.full<-prediction(test.probs,test$outcome)
pred.full<-performance(pred.full,"tpr","fpr")
plot(pred.full,main="ROC",col='darkblue',lwd=1.5)
#miRNA
test.probs=predict(mod2,type="response",newdata=test)
pred.full<-prediction(test.probs,test$outcome)
pred.full<-performance(pred.full,"tpr","fpr")
plot(pred.full,main="ROC",col='darkred',add=TRUE,lwd=1.5)
#Clinical
test.probs=predict(mod3,type="response",newdata=test)
pred.full<-prediction(test.probs,test$outcome)
pred.full<-performance(pred.full,"tpr","fpr")
plot(pred.full,main="ROC",col='darkgreen',add=TRUE,lwd=1.5)
legend(0.44,0.4,c("Nomogram(AUC=0.72, 0.68-0.77)","miRNA panel(AUC=0.64, 0.59-0.69)",
"Clinical(AUC=0.69, 0.64-0.73)"),col=c('darkblue','darkred','darkgreen'),lwd = 2,
cex = 0.5)
#trainging datasets
actualstest<-ifelse(spline$outcome=="1",1,0)
predtest<-predict(mod1,type="response")
predtest1<-predict(mod2,type="response")
predtest2<-predict(mod3,type="response")
roc1=roc(spline$outcome, predtest)#nomogram
ci(roc1,method="bootstrap")
roc2=roc(spline$outcome, predtest1)#miRNA
ci(roc2,method="bootstrap")
roc.test(roc1, roc2, method="bootstrap")
roc3=roc(spline$outcome, predtest2)#clinical
ci(roc3,method="bootstrap")
roc.test(roc1, roc3, method="bootstrap")
############################AUROC-training datasets
#nomogram
test.probs=predict(mod1,type="response")
pred.full<-prediction(test.probs,spline$outcome)
pred.full<-performance(pred.full,"tpr","fpr")
plot(pred.full,main="ROC",col='darkblue',lwd=1.5)
#miRNA
test.probs=predict(mod2,type="response")
pred.full<-prediction(test.probs,spline$outcome)
pred.full<-performance(pred.full,"tpr","fpr")
plot(pred.full,main="ROC",col='darkred',add=TRUE,lwd=1.5)
#Clinical
test.probs=predict(mod3,type="response")
pred.full<-prediction(test.probs,spline$outcome)
pred.full<-performance(pred.full,"tpr","fpr")
plot(pred.full,main="ROC",col='darkgreen',add=TRUE,lwd=1.5)
legend(0.44,0.4,c("Nomogram(AUC=0.74, 0.68-0.80)","miRNA panel(AUC=0.64, 0.57-0.71)",
                  "Clinical(AUC=0.68, 0.61-0.75)"),col=c('darkblue','darkred','darkgreen'),lwd = 2,
       cex = 0.5)
legend(0.44,0.4,c("Nomogram(AUC=0.70, 0.65-0.74)","miRNA panel(AUC=0.62, 0.57-0.67)",
  "Clinical(AUC=0.64, 0.59-0.69)"),col=c('darkblue','darkred','darkgreen')
  ,lwd = 2,cex = 0.5)
#testing datasets
#actualstrain<-ifelse(spline$outcome=="1",1,0)
actualstest<-ifelse(test$outcome=="1",1,0)
#predtrain<-predict(mod1,data=spline,type="response")
predtest<-predict(mod1,newdata=test,type="response")
predtest1<-predict(mod2,newdata=test,type="response")
predtest2<-predict(mod3,newdata=test,type="response")
#misClassError(actualstrain,predtrain)
library(plyr)
roc1=roc(test$outcome, predtest)#nomogram
ci(roc1,method="bootstrap")
roc2=roc(test$outcome, predtest1)#miRNA
ci(roc2,method="bootstrap")
roc.test(roc1, roc2, method="bootstrap")
#Vs clinical model
roc3=roc(test$outcome, predtest2)#clinical
ci(roc3,method="bootstrap")
roc.test(roc1, roc3, method="bootstrap")
#plotROC(actualstrain,predtrain)
plotROC(actualstest,predtest1)
###################################DCA curve
# add predicted values from model to data set
spline <-
  spline %>%
  mutate(
    Nomogram =
      broom::augment(mod1, type.predict = "response") %>%
      pull(".fitted")
  )
spline <-
  spline %>%
  mutate(
    miRNA.panels =
      broom::augment(mod2, type.predict = "response") %>%
      pull(".fitted")
  )
spline <-
  spline %>%
  mutate(
    Clinical =
      broom::augment(mod3, type.predict = "response") %>%
      pull(".fitted")
  )


#################################
setwd("//Users//houxiangqing//Desktop//temp//Figures//")
pdf("DCA_test.pdf",width = 8,height = 10)
dcurves::dca(outcome~Nomogram + miRNA.panels+Clinical,spline) %>% 
  plot(smooth=TRUE)
dev.off()
##########NR calculate
dcurves::dca(outcome ~ Nomogram,
    data = spline,
    as_probability = "Nomogram",
    
    label = list(Nomogram = "Nomogram")
) %>%
  net_intervention_avoided() %>%
  plot(smooth = TRUE)
