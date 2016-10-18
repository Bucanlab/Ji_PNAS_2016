#HIS & GHIS
dat=read.table(".../Essential_genes_master_table_IMPC_all_03282016_percFixed.txt",header=T,sep="\t",quote="",na.strings=c('.',"NA"))
datEG=dat[which(dat$EG_cell=="Y"),]
datNLG=dat[which(dat$NLG=="Y"),]

plotDat0=data.frame(dat$HIS_perc,rep('All',length(dat$HIS_perc)))
colnames(plotDat0) <- c("V1","V2")
plotDat1=data.frame(datEG$HIS_perc,rep('EG',length(datEG$HIS_perc)))
colnames(plotDat1) <- c("V1","V2")
plotDat2=data.frame(datNLG$HIS_perc,rep('NEG',length(datNLG$HIS_perc)))
colnames(plotDat2) <- c("V1","V2")


#plotDat=rbind(plotDat0,plotDat1,plotDat2,plotDat5,plotDat6)
plotDat=rbind(plotDat0,plotDat1,plotDat2)
colnames(plotDat) <- c("V1","Density")
plotDat$Density=factor(plotDat$Density,levels=c("EG","NEG","All"))
#plotDat$Density=factor(plotDat$Density,levels=c("EG","EG_PERI/POST","All","NEG","EG_PRE"))

test1=wilcox.test(plotDat1$V1,plotDat2$V1)
test1$p.value #2.454584e-122

plotDat0=data.frame(dat$GHIS_perc,rep('All',length(dat$GHIS_perc)))
colnames(plotDat0) <- c("V1","V2")
plotDat1=data.frame(datEG$GHIS_perc,rep('EG',length(datEG$GHIS_perc)))
colnames(plotDat1) <- c("V1","V2")
plotDat2=data.frame(datNLG$GHIS_perc,rep('NEG',length(datNLG$GHIS_perc)))
colnames(plotDat2) <- c("V1","V2")


#plotDat=rbind(plotDat0,plotDat1,plotDat2,plotDat5,plotDat6)
plotDat=rbind(plotDat0,plotDat1,plotDat2)
colnames(plotDat) <- c("V1","Density")
plotDat$Density=factor(plotDat$Density,levels=c("EG","NEG","All"))
#plotDat$Density=factor(plotDat$Density,levels=c("EG","EG_PERI/POST","All","NEG","EG_PRE"))

test1=wilcox.test(plotDat1$V1,plotDat2$V1)
test1$p.value #2.070937e-131


svg(".../figure1bc.svg",height=6,width=6)
ggplot(plotDat,aes(x=V1,group=Density,fill=as.factor(Density))) + 
  geom_histogram(aes(y = ..density..),position="identity", alpha=0.5,right=TRUE, binwidth=1) + 
  scale_fill_discrete(name="histogram") +
  geom_line(aes(x=V1,group=Density,color=Density), stat = 'density',position = "identity",adjust=0.5) +
  #scale_fill_manual(name="",values=c("#F8766D","#E76BF3","#E58700","#619CFF","##00BA38","black")) +
  theme_bw() + xlab("HIS percentile") + ylab("Density")+ xlim(c(0,100)) + ggtitle("") + ylim(c(0,0.03)) +
  theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=18),legend.position = c(0.85, 0.725))
dev.off()
