setwd("//Users//houxiangqing//Desktop//temp//Datasets_PPMI//")
install.packages("DecisionCurve")
install.packages("rms")
library(car)

library(pROC)
library(DecisionCurve)
spline <- read.csv("nomogram.csv")
spline$Education_Level<-factor(spline$Education_Level,labels = c('<12y','12-16y','>16y','Unknown'))
spline$Gender<-factor(spline$Gender,labels = c("Man","Woman"))
ddist <- datadist(spline)
options(datadist = "ddist")
fit <- lrm(outcome~Gender+Age+Education_Level+Transcriptional_Score, data=test, x=T, y=T)
summary(fit)

nom<- nomogram(fit, fun=plogis, fun.at=c(0.0001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.9999), lp=F, funlabel="Risk of PD") #构建Nomogram图
pdf("nomogram.pdf",width = 12,height = 8)
plot(nom)
dev.off()
/*Calibration curve*/
cal <- calibrate(fit, method="boot", B=1000)
plot(cal,
     xlim = c(0,1),
     xlab = "Nomogram-prediced PD risk (%)",
     ylab="Observed PD risk (%)",
     legend =FALSE,
     subtitles = FALSE)
abline(0,1,col="black",lty=2,lwd=2)
lines(cal[,c("predy","calibrated.orig")],type="l",lwd=2,col="red",pch=16)
lines(cal[,c("predy","calibrated.corrected")],type="l",lwd=2,col="green",pch=16)
legend(0.6,0.75,
       c("Ideal","Apparent","Bias-corrected"),
       lty = c(2,1,1),
       lwd = c(2,1,1),
       col = c("black","red","green"),
       bty="n") #"o"为加边框

cal<-calibrate(fit)
#plot(cal)
plot(cal,xlim=c(0.1,0.9),ylim=c(0.1,0.9),xlab = "Nomogram Predicted Risk", ylab = "Actual Risk")
####################Validation datasets
validate <- read.csv("validate.csv")
library(rms)
validate$Education_Level<-factor(validate$Education_Level,labels = c('<12y','12-16y','>16y','Unknown'))
validate$Gender<-factor(validate$Gender,labels = c("Man","Woman"))
ddist <- datadist(validate)
options(datadist = "ddist")

pred_f_validation<-predict(mylog,validate)
#(form<-as.formula(paste(names(data[5],"~",
#                             paste0(names(data)[1:4],collapse = "+")))))

fit.val <- lrm(outcome~Gender+Age+Education_Level+hsa_miR_4301+
                 hsa_miR_190a_5p+hsa_miR_22_3p+hsa_miR_3200_5p+
                 hsa_miR_3613_5p+hsa_miR_423_5p+hsa_miR_4433b_5p+
               +hsa_miR_4677_5p++hsa_miR_548b_5p+hsa_miR_654_5p, data=validate, x=T, y=T)

#pdf(file=paste(output_dir, "\\calibrate_testing.pdf", sep = ""),width=10,height=10) 
cal <- calibrate(fit.val,method="boot",B=1000)
plot(cal,xlab = "Nomogram-prediced PD risk (%)",ylab = "Observed PD risk (%)",)

############################train DCA
mod1 <- glm(outcome ~ Gender+Age+Education_Level+Transcriptional_Score, data = spline, family = binomial)
#miRNA model
mod2 <- glm(outcome ~ Transcriptional_Score, data = spline, family = binomial)
#Clinical model
mod3 <- glm(outcome ~ Gender+Age+Education_Level, data = spline, family = binomial)
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
pdf("DCA_train.pdf",width = 10,height = 10)
dcurves::dca(outcome~Nomogram + miRNA.panels+Clinical,spline) %>% 
  plot(smooth=TRUE)
dev.off()


############################validate DCA
mod1 <- glm(outcome~Gender+Age+Education_Level+hsa_miR_4301+
              hsa_miR_190a_5p+hsa_miR_22_3p+hsa_miR_3200_5p+
              hsa_miR_3613_5p+hsa_miR_423_5p+hsa_miR_4433b_5p+
              +hsa_miR_4677_5p++hsa_miR_548b_5p+hsa_miR_654_5p, data = spline, family = binomial)
#miRNA model
mod2 <- glm(outcome ~ hsa_miR_4301+
              hsa_miR_190a_5p+hsa_miR_22_3p+hsa_miR_3200_5p+
              hsa_miR_3613_5p+hsa_miR_423_5p+hsa_miR_4433b_5p+
              +hsa_miR_4677_5p++hsa_miR_548b_5p+hsa_miR_654_5p, data = spline, family = binomial)
#Clinical model
mod3 <- glm(outcome ~ Gender+Age+Education_Level, data = spline, family = binomial)
###################################DCA curve
# add predicted values from model to data set
validate <-
  validate %>%
  mutate(
    Nomogram =
      broom::augment(mod1, type.predict = "response") %>%
      pull(".fitted")
  )
validate <-
  validate %>%
  mutate(
    miRNA.panels =
      broom::augment(mod2, type.predict = "response") %>%
      pull(".fitted")
  )
validate <-
  validate %>%
  mutate(
    Clinical =
      broom::augment(mod3, type.predict = "response") %>%
      pull(".fitted")
  )


#################################
setwd("//Users//houxiangqing//Desktop//temp//Figures//")
pdf("DCA_validate.pdf",width = 10,height = 10)
dcurves::dca(outcome~Nomogram + miRNA.panels+Clinical,validate) %>% 
  plot(smooth=TRUE)
dev.off()



####################?ⲿ??֤??У׼????
external <- read.csv("D:\\desk\\tmp2.csv")
pred_f_external<-predict(mylog,external)

fit.vad<-lrm(external$group_spt~pred_f_external,
             data=external,x=T,y=T)

#pdf(file=paste(output_dir, "\\calibrate_testing.pdf", sep = ""),width=10,height=10) 
cal <- calibrate(fit.vad)
plot(cal)
#dev.off()
#####################################ROC
install.packages("ROCR")
library(ROCR)
spline <- read.csv("external1.csv")
spline1<-datadist(spline)
options(datadist='spline1')
bad.fit<-glm(outcome~hsa_miR_4301+hsa_miR_22_3p+hsa_miR_3200_5p+hsa_miR_3613_5p+hsa_miR_4433b_5p+hsa_miR_4677_5p+hsa_miR_654_5p+sex1+age+as.factor(Ethnicity1)+as.factor(Education_Level1),family = binomial,
             data=spline)
test.probs=predict(bad.fit,type = "response")
pred.full<-prediction(test.probs,spline$outcome)
pred.full<-performance(pred.full,"tpr","fpr")
plot(pred.full,main="ROC",col=1)

  
spline <- read.csv("external2.csv")
spline1<-datadist(spline)
options(datadist='spline1')
bad.fit<-glm(outcome~hsa_miR_4301+hsa_miR_22_3p+hsa_miR_3200_5p+hsa_miR_3613_5p+hsa_miR_4433b_5p+hsa_miR_4677_5p+hsa_miR_654_5p+sex1+age+as.factor(Ethnicity1)+as.factor(Education_Level1),family = binomial,
             data=spline)
test.probs=predict(bad.fit,type = "response")
pred.full<-prediction(test.probs,spline$outcome)
pred.full<-performance(pred.full,"tpr","fpr")
plot(pred.full,main="ROC",col=2,add=TRUE)

spline <- read.csv("external3.csv")
spline1<-datadist(spline)
options(datadist='spline1')
bad.fit<-glm(outcome~hsa_miR_4301+hsa_miR_22_3p+hsa_miR_3200_5p+hsa_miR_3613_5p+hsa_miR_4433b_5p+hsa_miR_4677_5p+hsa_miR_654_5p+sex1+age+as.factor(Ethnicity1)+as.factor(Education_Level1),family = binomial,
             data=spline)
test.probs=predict(bad.fit,type = "response")
pred.full<-prediction(test.probs,spline$outcome)
pred.full<-performance(pred.full,"tpr","fpr")
plot(pred.full,main="ROC",col=3,add=TRUE)

spline <- read.csv("external4.csv")
spline1<-datadist(spline)
options(datadist='spline1')
bad.fit<-glm(outcome~hsa_miR_4301+hsa_miR_22_3p+hsa_miR_3200_5p+hsa_miR_3613_5p+hsa_miR_4433b_5p+hsa_miR_4677_5p+hsa_miR_654_5p+sex1+age+as.factor(Ethnicity1)+as.factor(Education_Level1),family = binomial,
             data=spline)
test.probs=predict(bad.fit,type = "response")
pred.full<-prediction(test.probs,spline$outcome)
pred.full<-performance(pred.full,"tpr","fpr")
plot(pred.full,main="ROC",col=4,add=TRUE)

spline <- read.csv("external5.csv")
spline1<-datadist(spline)
options(datadist='spline1')
bad.fit<-glm(outcome~hsa_miR_4301+hsa_miR_22_3p+hsa_miR_3200_5p+hsa_miR_3613_5p+hsa_miR_4433b_5p+hsa_miR_4677_5p+hsa_miR_654_5p+sex1+age+as.factor(Ethnicity1)+as.factor(Education_Level1),family = binomial,
             data=spline)
test.probs=predict(bad.fit,type = "response")
pred.full<-prediction(test.probs,spline$outcome)
pred.full<-performance(pred.full,"tpr","fpr")
plot(pred.full,main="ROC",col=5,add=TRUE)
legend(c("BL(AUC=0.75)","M6(AUC=0.81)","M12(AUC=0.72)",
                   "M24(AUC=0.70)","M36(AUC=0.69)"))


spline <- read.csv("D:\\Desk\\work\\work41_luowt\\paper1_predict IgE mediated allergy\\Figure1\\tmp0823\\pre_figure\\external2.csv")
spline1<-datadist(spline)
options(datadist='spline1')
bad.fit<-glm(group_spt~sex+minzu+Q6+q40_grp+Q47+C1+C2+C4+C11,family = binomial,
             data=spline)
test.probs=predict(bad.fit,type = "response")
pred.full<-prediction(test.probs,spline$group_spt)
pred.full<-performance(pred.full,"tpr","fpr")
plot(pred.full,main="ROC",col=3,add=TRUE)
#####################DCA
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
mod2 <- glm(outcome ~ Gender + Age + Education_Level+Transcriptional_Score, data = spline, family = binomial)
#Summary tbl_regression(mod2, exponentiate = TRUE)
mod3 <- glm(outcome ~ Transcriptional_Score, data = spline, family = binomial)
tbl_regression(mod2, exponentiate = TRUE)
spline <-
  spline %>%
  mutate(
    nomogram =
      broom::augment(mod2, type.predict = "response") %>%
      pull(".fitted")
  )
spline <-
  spline %>%
  mutate(
    miRNA =
      broom::augment(mod3, type.predict = "response") %>%
      pull(".fitted")
  )
dca(outcome ~ miRNA + nomogram,
    data = spline,
    thresholds = seq(0, 1.0, 0.1)
) %>%
  plot(smooth = TRUE)

dca(outcome ~ miRNA + nomogram,
    data = spline,
    thresholds = seq(0, 1, 0.1),
    #label = list(c(miRNA="miRNA panels",nomogram = "Nomogram"))
) %>%
  plot(smooth = TRUE)
library(DCA)
model <- predict(lrm(outcome ~ nomogram), type="fitted") 

fit <- lrm(outcome~Gender+Age+Education_Level+Transcriptional_Score, data=spline, x=T, y=T)

fit1 <- lrm(outcome~Transcriptional_Score, data=spline, x=T, y=T)

fit2 <- lrm(outcome~Gender+Age+Education_Level, data=spline, x=T, y=T)

library(dcurves)
library(gtsummary)
library(dplyr)
library(tidyr)
dcurves::dca(outcome~nomogram+miRNA,spline) %>% plot(smooth=TRUE)
