library(ggplot2)
library(reshape2)

# SFARI score
setwd("...")
dat=read.table("SFARI_table_all.txt",header=T,sep="\t")
dat1=dat[c(1:7),c(1:4)]

dat2=melt(dat1,id.vars='SFARI_score')
colnames(dat2)[2]="GeneSet"
levels(dat2$SFARI_score)=c("1+1S\n(n=16)","2+2S\n(n=37)","3+3S\n(n=109)","4+4S\n(n=208)","5\n(n=115)","6\n(n=19)","S\n(n=39)")
dat2$SFARI_score=factor(dat2$SFARI_score,levels(dat2$SFARI_score)[c(7,1:6)])
print(levels(dat2$SFARI_score))

png("barplot_SFARIscore.png",height=6,width=7,units='in',res=500)
ggplot(dat2,aes(SFARI_score,value,fill=GeneSet))+
	scale_fill_manual(values = c("#F8766D","#00BFC4", "grey"))+
	geom_bar(position = "dodge", stat="identity")+
	theme_bw() + xlab("SFARI Score Category") + ylab("Proportion")+ ggtitle("") + scale_y_continuous(limits = c(0, 0.8))+
	theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=13),axis.title=element_text(size=18))
dev.off()

#ClinGen HI
setwd("...")
dat=read.table("clinGen_HI_table.txt",header=T,sep="\t")
dat1=dat[c(1:5),c(1:4)]

dat2=melt(dat1,id.vars='ClinGen_HI_score')
colnames(dat2)[2]="GeneSet"
dat2$ClinGen_HI_score = factor(dat2$ClinGen_HI_score,levels(dat2$ClinGen_HI_score)[c(4,3,2,1,5)])
levels(dat2$ClinGen_HI_score)=c("Sufficient evidence\n(n=239)","Some evidence\n(n=41)","Little evidence\n(n=47)","No evidence\n(n=200)","Recessive/Not HI\n(n=89)")
print(levels(dat2$ClinGen_HI_score))

png("barplot_ClinGen_HI_score.png",height=6,width=7,units='in',res=500)
ggplot(dat2,aes(ClinGen_HI_score,value,fill=GeneSet))+
	scale_fill_manual(values = c("#F8766D","#00BFC4", "grey"))+
	geom_bar(position = "dodge", stat="identity")+
	theme_bw() + xlab("ClinGen Haploinsufficiency Rating") + ylab("Proportion")+ ggtitle("") + scale_y_continuous(limits = c(0, 0.8))+
	theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=13),axis.title=element_text(size=18))
dev.off()
