library(ggplot2)
dat=read.table(".../masterPheno_indiv.txt",header=T,sep="\t")

plotDat1=data.frame(dat$Proband.Verbal.IQ,dat$Proband.SRS.Total)
plotDat2=data.frame(dat$Proband.Nonverbal.IQ,dat$Proband.SRS.Total)
plotDat3=data.frame(dat$Proband.Full.Scale.IQ,dat$Proband.SRS.Total)

cor(plotDat1[,1],plotDat1[,2],use="complete.obs")
png(".../xyPlot_SRS_verbalIQ.png",height=6,width=6,units='in',res=500)
ggplot(plotDat1,aes(x=dat.Proband.Verbal.IQ,y=dat.Proband.SRS.Total)) + 
  geom_point(shape=1) +
  theme_bw() + xlab("Proband Verbal IQ") + ylab("Proband SRS Total") + ggtitle("") +
  theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=13),axis.title=element_text(size=18))
dev.off()

cor(plotDat2[,1],plotDat2[,2],use="complete.obs")
png(".../xyPlot_SRS_nonVerbalIQ.png",height=6,width=6,units='in',res=500)
ggplot(plotDat2,aes(x=dat.Proband.Nonverbal.IQ,y=dat.Proband.SRS.Total)) + 
  geom_point(shape=1) +
  theme_bw() + xlab("Proband Nonverbal IQ") + ylab("Proband SRS Total") + ggtitle("") +
  theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=13),axis.title=element_text(size=18))
dev.off()

cor(plotDat1[,1],plotDat1[,2],use="complete.obs")
png(".../xyPlot_SRS_fullScaleIQ.png",height=6,width=6,units='in',res=500)
ggplot(plotDat1,aes(x=dat.Proband.Verbal.IQ,y=dat.Proband.SRS.Total)) + 
  geom_point(shape=1) +
  theme_bw() + xlab("Proband Verbal IQ") + ylab("Proband SRS Total") + ggtitle("") +
  theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=13),axis.title=element_text(size=18))
dev.off()

#plot only IQ<=50107
plotDat1_LF=plotDat1[which(plotDat1$dat.Proband.Verbal.IQ<=50),]
plotDat2_LF=plotDat2[which(plotDat2$dat.Proband.Nonverbal.IQ<=50),]
plotDat3_LF=plotDat3[which(plotDat3$dat.Proband.Full.Scale.IQ<=50),]

cor(plotDat1_LF[,1],plotDat1_LF[,2],use="complete.obs") #-0.2353728
png(".../xyPlot_SRS_verbalIQ_IQbelow50.png",height=6,width=6,units='in',res=500)
ggplot(plotDat1_LF,aes(x=dat.Proband.Verbal.IQ,y=dat.Proband.SRS.Total)) + 
  geom_point(shape=1) +
  theme_bw() + xlab("Proband Verbal IQ") + ylab("Proband SRS Total") + ggtitle("") +
  theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=13),axis.title=element_text(size=18))
dev.off()

cor(plotDat2_LF[,1],plotDat2_LF[,2],use="complete.obs") #-0.2133093
png(".../xyPlot_SRS_nonVerbalIQ_IQbelow50.png",height=6,width=6,units='in',res=500)
ggplot(plotDat2_LF,aes(x=dat.Proband.Nonverbal.IQ,y=dat.Proband.SRS.Total)) + 
  geom_point(shape=1) +
  theme_bw() + xlab("Proband Nonverbal IQ") + ylab("Proband SRS Total") + ggtitle("") +
  theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=13),axis.title=element_text(size=18))
dev.off()

cor(plotDat3_LF[,1],plotDat3_LF[,2],use="complete.obs") #-0.2447624
png(".../xyPlot_SRS_fullScaleIQ_IQbelow50.png",height=6,width=6,units='in',res=500)
ggplot(plotDat3_LF,aes(x=dat.Proband.Full.Scale.IQ,y=dat.Proband.SRS.Total)) + 
  geom_point(shape=1) +
  theme_bw() + xlab("Proband Full Scale IQ") + ylab("Proband SRS Total") + ggtitle("") +
  theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=13),a046xis.title=element_text(size=18))
dev.off()

#plot only IQ>50
plotDat1_HF=plotDat1[which(plotDat1$dat.Proband.Verbal.IQ>50),]
plotDat2_HF=plotDat2[which(plotDat2$dat.Proband.Nonverbal.IQ>50),]
plotDat3_HF=plotDat3[which(plotDat3$dat.Proband.Full.Scale.IQ>50),]

cor(plotDat1_HF[,1],plotDat1_HF[,2],use="complete.obs") #-0.04566907
png(".../xyPlot_SRS_verbalIQ_IQabove50.png",height=6,width=6,units='in',res=500)
ggplot(plotDat1_HF,aes(x=dat.Proband.Verbal.IQ,y=dat.Proband.SRS.Total)) + 
  geom_point(shape=1) +
  theme_bw() + xlab("Proband Verbal IQ") + ylab("Proband SRS Total") + ggtitle("") +
  theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=13),axis.title=element_text(size=18))
dev.off()

cor(plotDat2_HF[,1],plotDat2_HF[,2],use="complete.obs") #-0.1323556
png(".../xyPlot_SRS_nonVerbalIQ_IQabove50.png",height=6,width=6,units='in',res=500)
ggplot(plotDat2_HF,aes(x=dat.Proband.Nonverbal.IQ,y=dat.Proband.SRS.Total)) + 
  geom_point(shape=1) +
  theme_bw() + xlab("Proband Nonverbal IQ") + ylab("Proband SRS Total") + ggtitle("") +
  theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=13),axis.title=element_text(size=18))
dev.off()

cor(plotDat3_HF[,1],plotDat3_HF[,2],use="complete.obs") #-0.1069003
png(".../xyPlot_SRS_fullScaleVerbalIQ_IQabove50.png",height=6,width=6,units='in',res=500)
ggplot(plotDat3_HF,aes(x=dat.Proband.Full.Scale.IQ,y=dat.Proband.SRS.Total)) + 
  geom_point(shape=1) +
  theme_bw() + xlab("Proband Full Scale IQ") + ylab("Proband SRS Total") + ggtitle("") +
  theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=13),axis.title=element_text(size=18))
dev.off()
