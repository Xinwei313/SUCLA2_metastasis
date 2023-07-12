

###### Survival package
install.packages("lattice")
library(lattice)


######## plot KM curves
install.packages("reshape2")
library(reshape2)
install.packages("data.table")
library(data.table)
install.packages("zoo")
library("zoo")
install.packages("survminer")
library("survminer")
install.packages("survival")
library("survival")

#### ROC package
install.packages("gplots")
library(gplots)
install.packages("ROCR")
library(ROCR)
install.packages("dplyr")
library("dplyr")
install.packages("lubridate")
library(lubridate)
install.packages("robustbase")
install.packages("caret")
library(caret)
install.packages("pROC")
library(pROC)

install.packages("timeROC")
install.packages("prodlim")
library(prodlim)
library(quantreg)
install.packages("quantreg")
install.packages("polspline")
library("polspline")
library(timeROC)
install.packages("foreign")
library(foreign)

# install.packages("rms")
library(rms)


install.packages("glmnet")
library(glmnet)  # manually load this package by downloading and installing from Toll. 


#=====================================================================================
#
#  Code chunk 1   
#
#=====================================================================================


#This RData is directly downloaded from GEO with accession number GSE7390
load("GSE7390_transbig2006affy.RData")

par(mar=c(4,4,2,1))
title <- paste ("GSE7390")
plotDensities(data, main=title, legend=F)

Markergenes = c( "SUCLA2", "USP10")  
Annotation = annot[annot$NCBI.gene.symbol %in% Markergenes,]

Markers_expr = data[,c("202930_s_at",                                  # "SUCLA2"
                       "209136_s_at", "209137_s_at")]                  # "USP10"

USP10_expr = Markers_expr[,c(2:3)]
USP10_expr = rowSums(USP10_expr)/2


Markers_expr = cbind(Markers_expr[,1], USP10_expr, Markers_expr[,1]*USP10_expr)
colnames(Markers_expr) = c("SUCLA2", "USP10", "SUCLA2*USP10")


Clinical_info = demo[rownames(Markers_expr),c(5,6,12,14:19)]   

time = as.numeric(Clinical_info[,8])     #DMFS
status = as.numeric(Clinical_info[,9])


predicted = as.numeric(Markers_expr[,1])  # "SUCLA2"

predicted = as.numeric(Markers_expr[,2])  # "USP10" 

predicted = as.numeric(Markers_expr[,3])  # "SUCLA2*USP10"


predicted0=as.numeric(predicted)
predicted0=predicted[!is.na(status)]
status0=na.omit(status)
pred <- prediction(as.numeric(predicted0),status0)
perf <- performance(pred,"tpr","fpr") 
performance(pred,"auc") 
dev.new()
plot(perf,colorize=FALSE, col="blue") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc1=roc(status0,as.numeric(predicted0))
AUC=roc1$auc
AUC
## optimal combination
opt <- which.max(rowSums(cbind(roc1$sensitivities,roc1$specificities)))
## optimal cut-off point 
sort(predicted0,F)[opt]


groups=matrix(0,1,length(status0))
groups[predicted0<=sort(predicted0,F)[opt]]=1
groups[predicted0>sort(predicted0,F)[opt]]=2
groups=t(groups)
groups=as.numeric(groups)


fit<- survfit(Surv(time/365, status) ~ groups, data = as.data.frame(predicted))

ggsurv <- ggsurvplot(fit, pval = TRUE, conf.int = FALSE,
                     palette = c("#68BFCF","#F94040"),
                     xlab='Time (years)', 
                     ylab='DMFS probability',
                     ylim=c(0,1), break.y.by=0.5,
                     font.x=c(24,"bold"), font.y=c(24,"bold"), font.tickslab=c(24, "bold"),
                     font.main=c(18,"bold"),
                     risk.table = F, legend=c(0.4,0.2),
                     legend.title=c(""),
                     legend.labs=c("SUCLA2 Low (n = 94)","SUCLA2 High (n = 104)"),
                     #font.legend=c(20),
                     pval.coord=c(2,0.35), pval.size = 10)
dev.new()
ggsurv


fit<- survfit(Surv(time/365, status) ~ groups, data = as.data.frame(predicted))

ggsurv <- ggsurvplot(fit, pval = TRUE, conf.int = FALSE,
                     palette = c("#68BFCF","#F94040"),
                     xlab='Time (years)', 
                     ylab='DMFS probability',
                     ylim=c(0,1), break.y.by=0.5,
                     font.x=c(24,"bold"), font.y=c(24,"bold"), font.tickslab=c(24, "bold"),
                     font.main=c(18,"bold"),
                     risk.table = F, legend=c(0.4,0.2),
                     legend.title=c(""),
                     legend.labs=c("USP10 Low (n = 86)","USP10 High (n = 112)"),
                     #font.legend=c(20),
                     pval.coord=c(2,0.35), pval.size = 10)
dev.new()
ggsurv


fit<- survfit(Surv(time/365, status) ~ groups, data = as.data.frame(predicted))

ggsurv <- ggsurvplot(fit, pval = TRUE, conf.int = FALSE,
                     palette = c("#68BFCF","#F94040"),
                     xlab='Time (years)', 
                     ylab='DMFS probability',
                     ylim=c(0,1), break.y.by=0.5,
                     font.x=c(24,"bold"), font.y=c(24,"bold"), font.tickslab=c(24, "bold"),
                     font.main=c(18,"bold"),
                     risk.table = F, legend=c(0.4,0.2),
                     legend.title=c(""),
                     legend.labs=c("SUCLA2*USP10 Low (n = 117)","SUCLA2*USP10 High (n = 81)"),
                     #font.legend=c(20),
                     pval.coord=c(2,0.35), pval.size = 10)
dev.new()
ggsurv








