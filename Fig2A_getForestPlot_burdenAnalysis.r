#forest plot
setwd("...")
library(forestplot)
meanD=c(0.1169762,0.03003932,0.07940217,0.03743397,0.04199952,-0.01503358)
lowCL=c(0.07154123,-0.01574115,0.03238893,-0.01002884,-0.00410884,-0.06180866)
highCL=c(0.1595998,0.07744043,0.1273962,0.08266127,0.08866722,0.03079507)
row_names=c("dnLoF","dnNSD","inhRD")
meanD=data.frame(c(0.1169762,0.07940217,0.04199952),c(0.03003932,0.03743397,-0.01503358))
colnames(meanD)=c("EG","NEG")
lowCL=data.frame(c(0.07154123,0.03238893,-0.00410884),c(-0.01574115,-0.01002884,-0.06180866))
colnames(lowCL)=c("EG","NEG")
highCL=data.frame(c(0.1595998,0.1273962,0.08866722),c(0.07744043,0.08266127,0.03079507))
colnames(highCL)=c("EG","NEG")


svg("forestPlot_CohensD_largeFont.svg",height=6,width=6)
formatSetting=fpTxtGp()
formatSetting$xlab$cex=2.5
formatSetting$ticks$cex=2.0
formatSetting$label$cex=2.5
formatSetting$legend$cex=2.5

forestplot(row_names, meanD, lowCL, highCL, is.summary=F ,col=fpColors(box=c("#F8766D","#00BFC4"),line=c("#F8766D","#00BFC4"),summary=c("#F8766D","#00BFC4")),clip=c(-2,2),xlab="Effect size",ci.vertices=T,ci.vertices.height=0.18 ,txt_gp=formatSetting,lty.ci = 1, lwd.ci=2.5, boxsize=0.30, legend=c("EG","NEG"))
dev.off()
