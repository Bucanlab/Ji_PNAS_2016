library(ggplot2)

datEG=read.table(".../Sanders_2015_S6_TADA.txt.EG_all_cell_v2_100000",header=T,sep='\t')
datNEG=read.table(".../Sanders_2015_S6_TADA.txt.NEG_all_v2_100000",header=T,sep='\t')
datAll=read.table(".../Sanders_2015_S6_TADA.txt.allGenes_100000",header=T,sep='\t')


#log(-log(q))
plotDat1=data.frame(log10(-log10(datEG$TADA_FDR)),log10(-log10(datEG$TADA_FDR_expected)),rep('EG',length(datEG$TADA_FDR)))
colnames(plotDat1) <- c("V1","V2")
plotDat2=data.frame(log10(-log10(datNEG$TADA_FDR)),log10(-log10(datNEG$TADA_FDR_expected)),rep('NEG',length(datNEG$TADA_FDR)))
colnames(plotDat2) <- c("V1","V2")
plotDat3=data.frame(log10(-log10(datAll$TADA_FDR_expected_2.5)),log10(-log10(datAll$TADA_FDR_expected)),rep('2.5perc',length(datAll$TADA_FDR)))
colnames(plotDat3) <- c("V1","V2")
plotDat4=data.frame(log10(-log10(datAll$TADA_FDR_expected_97.5)),log10(-log10(datAll$TADA_FDR_expected)),rep('97.5perc',length(datAll$TADA_FDR)))
colnames(plotDat4) <- c("V1","V2")

plotDat_main=rbind(plotDat1,plotDat2)
colnames(plotDat_main) <- c("TADA_FDR","TADA_FDR_expected","GeneSet")
plotDat_main$GeneSet=factor(plotDat_main$GeneSet,levels=c("EG","NEG"))

plotDat_cb=rbind(plotDat3,plotDat4)
colnames(plotDat_cb) <- c("TADA_FDR","TADA_FDR_expected","GeneSet")
plotDat_cb$GeneSet=factor(plotDat_cb$GeneSet,levels=c("2.5perc","97.5perc"))


png(".../qqPlot_TADA_EGcell_NEG_v2_confidenceBand_perm100000.png",height=6,width=6,units='in',res=500)
ggplot() + 
  #geom_line(data=plotDat_cb, aes(x=TADA_FDR_expected,y=TADA_FDR,group=GeneSet), linetype="dotted" ) +
  geom_ribbon(data=datAll,aes(x=log10(-log10(TADA_FDR_expected)), ,ymin=log10(-log10(TADA_FDR_expected_2.5)), ymax=log10(-log10(TADA_FDR_expected_97.5))),alpha=0.2) +
  geom_abline(slope=1,linetype="dashed") +
  geom_abline(intercept=0,slope=0, linetype="dashed",color="red") +
  geom_abline(intercept=-0.5213902,slope=0, linetype="dashed",color="blue") +
  geom_point(data=plotDat_main, aes(x=TADA_FDR_expected,y=TADA_FDR,group=GeneSet,color=GeneSet),shape=16,size=2.3,alpha=0.8) + 
  scale_fill_discrete(name="") +
  theme_bw() + xlab("Expected TADA FDR per gene") + ylab("Observed TADA FDR per gene") + ggtitle("") +# xlim(c(NA,1))+ ylim(c(NA,1)) +labels=c("0.1","1","10") 
  scale_x_continuous(limits=c(NA,1),labels=c("","0.8","0.1","1E-10")) +
  scale_y_continuous(limits=c(NA,1),labels=c("","0.8","0.1","1E-10")) +
  theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=13),axis.title=element_text(size=18),legend.position=c(0.87,0.18))
dev.off()

#Together

datEGNEG=read.table(".../Sanders_2015_S6_TADA.txt.EG_NEG_all_100000",header=T,sep='\t')
plotDat=data.frame(log10(-log10(datEGNEG$TADA_FDR)),log10(-log10(datEGNEG$TADA_FDR_expected)),datEGNEG$glName)
colnames(plotDat) <- c("TADA_FDR","TADA_FDR_expected","GeneSet")

png(".../qqPlot_TADA_EG_NEG_perm100000.png",height=9,width=9,units='in',res=500)
ggplot(plotDat, aes(x=TADA_FDR_expected,y=TADA_FDR,group=GeneSet,fill=GeneSet)) + 
  geom_point(shape=21,size=6) + 
  geom_abline(slope=1,linetype="dashed") +
  geom_abline(intercept=0,slope=0, linetype="dashed",color="red") +
  scale_fill_discrete(name="") +
  theme_bw() + xlab("Expected TADA -log(q) per gene") + ylab("Observed TADA -log(q) per gene") + ggtitle("") +# xlim(c(NA,1))+ ylim(c(NA,1)) +labels=c("0.1","1","10") 
  scale_x_continuous(limits=c(NA,1),labels=c("","0.1","1","10")) +
  scale_y_continuous(limits=c(NA,1),labels=c("","0.1","1","10")) +
  theme(plot.title=element_text(size=20),legend.text=element_text(size=18),axis.text=element_text(size=18),axis.title=element_text(size=18),legend.position=c(0.9,0.1))
dev.off()
