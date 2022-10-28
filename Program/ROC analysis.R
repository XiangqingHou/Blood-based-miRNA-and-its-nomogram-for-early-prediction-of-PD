
################################import data
setwd("//Users//houxiangqing//Desktop//temp//Datasets_PPMI//")
spline <- read.csv("train.csv")
validation <- read.csv("validation.csv")
test <- read.csv("External1.csv")
test$Education_Level<-factor(test$Education_Level,labels = c('<12y','12-16y','>16y','Unknown'))
test$Gender<-factor(test$Gender,labels = c("Man","Woman"))
test1<-datadist(test)
options(datadist='test1')
spline$Education_Level<-factor(spline$Education_Level,labels = c('<12y','12-16y','>16y','Unknown'))
spline$Gender<-factor(spline$Gender,labels = c("Man","Woman"))
spline1<-datadist(spline)
options(datadist='spline1')
validation$Education_Level<-factor(validation$Education_Level,labels = c('<12y','12-16y','>16y','Unknown'))
validation$Gender<-factor(validation$Gender,labels = c("Man","Woman"))
validation1<-datadist(validation)
options(datadist='validation1')

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
#legend of training sets
legend(0.44,0.4,c("Nomogram(AUC=0.70, 0.65-0.74)","miRNA panel(AUC=0.62, 0.57-0.67)",
"Clinical(AUC=0.64, 0.59-0.69)"),col=c('darkblue','darkred','darkgreen')
,lwd = 2,cex = 0.5)
############################AUROC-validation datasets
#nomogram
test.probs=predict(mod1,type="response")
pred.full<-prediction(test.probs,validation$outcome)
pred.full<-performance(pred.full,"tpr","fpr")
plot(pred.full,main="ROC",col='darkblue',lwd=1.5)
#miRNA
test.probs=predict(mod2,type="response")
pred.full<-prediction(test.probs,validation$outcome)
pred.full<-performance(pred.full,"tpr","fpr")
plot(pred.full,main="ROC",col='darkred',add=TRUE,lwd=1.5)
#Clinical
test.probs=predict(mod3,type="response")
pred.full<-prediction(test.probs,validation$outcome)
pred.full<-performance(pred.full,"tpr","fpr")
plot(pred.full,main="ROC",col='darkgreen',add=TRUE,lwd=1.5)
#legend of validation sets
legend(0.44,0.4,c("Nomogram(AUC=0.74, 0.68-0.80)","miRNA panel(AUC=0.64, 0.57-0.71)",
                  "Clinical(AUC=0.68, 0.61-0.75)"),col=c('darkblue','darkred','darkgreen'),lwd = 2,
       cex = 0.5)