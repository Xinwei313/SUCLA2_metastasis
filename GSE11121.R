

##### GSE12093 #####


library(GEOquery)
library(limma)
library(umap)


#=====================================================================================
#
#  Code chunk 1 Normalization
#
#=====================================================================================


GSE11121 <- getGEO('GSE11121', GSEMatrix = TRUE, getGPL = FALSE)
exprSet <- exprs(GSE11121[[1]])

# log2 transform
exprSet <- log2(exprSet)

par(mar=c(4,4,2,1))
title <- paste ("GSE11121", "/", annotation(GSE11121[[1]]), " value distribution", sep ="")
plotDensities(exprSet, main=title, legend=F)


#=====================================================================================
#
#  Code chunk 2 Annotation
#
#=====================================================================================


#Annotation
gpl96 <- getGEO('GPL96', destdir = '.')

Annotation = Table(gpl96)

Markergenes = c( "SUCLA2", "USP10")  

Annotation = Annotation[Annotation$`Gene Symbol` %in% Markergenes,]
Annotation = Annotation[,c(1,11)]

exprSet = exprSet[rownames(exprSet) %in% Annotation$ID,]


#Take average of the repeated ID
SUCLA2_expr = exprSet[rownames(exprSet) %in% Annotation$ID[Annotation$`Gene Symbol`=="SUCLA2"],]
#SUCLA2_expr = colSums(SUCLA2_expr)/nrow(SUCLA2_expr)

USP10_expr = exprSet[rownames(exprSet) %in% Annotation$ID[Annotation$`Gene Symbol`=="USP10"],]
USP10_expr = colSums(USP10_expr)/nrow(USP10_expr)

Markers_expr = rbind(SUCLA2_expr, USP10_expr, SUCLA2_expr*USP10_expr)
rownames(Markers_expr) = c("SUCLA2", "USP10", "SUCLA2*USP10")


#=====================================================================================
#
#  Code chunk 3 Clinical_info
#
#=====================================================================================


Clinical_info = pData(GSE11121[[1]])
colnames(Clinical_info)

Clinical_info = Clinical_info[colnames(Markers_expr),c(14,15,44,49)]   

time = as.numeric(Clinical_info[,4])     
status = as.numeric(Clinical_info[,3])


#=====================================================================================
#
#  Code chunk 4 K-M analysis
#
#=====================================================================================


predicted = as.numeric(Markers_expr[1,])  # "SUCLA2"

predicted = as.numeric(Markers_expr[2,])  # "USP10" 

predicted = as.numeric(Markers_expr[3,])  # "SUCLA2*USP10"


predicted0=as.numeric(predicted)
predicted0=predicted[!is.na(status)]
status0=na.omit(status)
pred <- prediction(as.numeric(predicted0),status0)
perf <- performance(pred,"tpr","fpr")  
performance(pred,"auc") 
dev.new()
plot(perf,colorize=FALSE, col="blue") # plot ROC curve
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


fit<- survfit(Surv(time/12, status) ~ groups, data = as.data.frame(predicted))

ggsurv <- ggsurvplot(fit, pval = TRUE, conf.int = FALSE,
                     palette = c("#68BFCF","#F94040"),
                     xlab='Time (years)', 
                     ylab='DMFS probability',
                     ylim=c(0,1), break.y.by=0.5,
                     font.x=c(24,"bold"), font.y=c(24,"bold"), font.tickslab=c(24, "bold"),
                     font.main=c(18,"bold"),
                     risk.table = F, legend=c(0.4,0.2),
                     legend.title=c(""),
                     legend.labs=c("SUCLA2 Low (n = 169)","SUCLA2 High (n = 31)"),
                     #font.legend=c(20),
                     pval.coord=c(2,0.35), pval.size = 10)
dev.new()
ggsurv


fit<- survfit(Surv(time/12, status) ~ groups, data = as.data.frame(predicted))

ggsurv <- ggsurvplot(fit, pval = TRUE, conf.int = FALSE,
                     palette = c("#68BFCF","#F94040"),
                     xlab='Time (years)', 
                     ylab='DMFS probability',
                     ylim=c(0,1), break.y.by=0.5,
                     font.x=c(24,"bold"), font.y=c(24,"bold"), font.tickslab=c(24, "bold"),
                     font.main=c(18,"bold"),
                     risk.table = F, legend=c(0.4,0.2),
                     legend.title=c(""),
                     legend.labs=c("USP10 Low (n = 72)","USP10 High (n = 128)"),
                     #font.legend=c(20),
                     pval.coord=c(2,0.35), pval.size = 10)
dev.new()
ggsurv


fit<- survfit(Surv(time/12, status) ~ groups, data = as.data.frame(predicted))

ggsurv <- ggsurvplot(fit, pval = TRUE, conf.int = FALSE,
                     palette = c("#68BFCF","#F94040"),
                     xlab='Time (years)', 
                     ylab='DMFS probability',
                     ylim=c(0,1), break.y.by=0.5,
                     font.x=c(24,"bold"), font.y=c(24,"bold"), font.tickslab=c(24, "bold"),
                     font.main=c(18,"bold"),
                     risk.table = F, legend=c(0.4,0.2),
                     legend.title=c(""),
                     legend.labs=c("SUCLA2*USP10 Low (n = 49)","SUCLA2*USP10 High (n = 151)"),
                     #font.legend=c(20),
                     pval.coord=c(2,0.35), pval.size = 10)
dev.new()
ggsurv

