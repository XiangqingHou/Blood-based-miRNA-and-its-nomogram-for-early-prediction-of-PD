setwd("//Users//houxiangqing//Desktop//temp//Datasets_PPMI//")
# install dcurves to perform DCA from CRAN
install.packages("dcurves")
install.packages("rms")
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
#Calibration curve
cal <- calibrate(fit, method="boot", B=1000)
plot(cal,xlim=c(0.1,0.9),ylim=c(0.1,0.9),xlab = "Nomogram Predicted Risk", ylab = "Actual Risk")
############################nomogram DCA
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
##########NR calculate
dcurves::dca(outcome ~ Nomogram,
    data = spline,
    as_probability = "Nomogram",
    
    label = list(Nomogram = "Nomogram")
) %>%
  net_intervention_avoided() %>%
  plot(smooth = TRUE)