library(ggplot2)
library(plyr)

setwd("...")

dat=read.table("enrichmentPerModule_EGcell_NEG_v2_TADA0.5_updated.txt",header=T,sep='\t')
dat_EG=dat[,c(1,2,3,4,5)]
dat_EG=data.frame(dat_EG,rep(NA,dim(dat_EG)[1]))
names(dat_EG)=c("Module","EG_all","NLG_all","OR","P","Category")
dat_EG$Category[which(dat_EG$OR<1)] = "3NEG enrichment"
dat_EG$Category[which(dat_EG$OR>1)] = "1EG enrichment"
FDR_EG=p.adjust(dat_EG$P,method="BH")

plotDat1=data.frame(substr(dat_EG$Module,1,3),-log10(dat_EG$P),-log10(FDR_EG),dat_EG$Category)
names(plotDat1)=c("Module","Pval","FDR","Category")

dat_TADA=dat[,c(1,6,7,8,9)]
dat_TADA=data.frame(dat_TADA,rep('2Potential ASD gene enrichment',dim(dat_EG)[1]))
names(dat_TADA)=c("Module","TADA0.5","nonTADA0.5","OR","P","Category")
FDR_TADA=p.adjust(dat_TADA$P,method="BH")

plotDat2=data.frame(substr(dat_TADA$Module,1,3),log10(dat_TADA$P),log10(FDR_TADA),dat_TADA$Category)
names(plotDat2)=c("Module","Pval","FDR","Category")

plotDat=rbind(plotDat1,plotDat2)
names(plotDat)=c("Module","Pval","FDR","Category")
plotDat$Category=factor(plotDat$Category,levels=c("1EG enrichment","2Potential ASD gene enrichment","3NEG enrichment"))
orderedEGmods=substr(c("M01","M02","M07","M26","M22","M08","M37","M17","M25","M34","M14","M16","M10","M15","M18","M24","M33"),1,3)
orderedNEGmods=substr(c("M04","M19","M12","M11","M06","M38","M09","M13","M30","M03","M32","M23","M29","M35","M27","M41","M05","M20","M36","M40","M31","M39","M21","M28"),1,3)
plotDat$Module=factor(plotDat$Module,levels=c(rev(orderedEGmods),orderedNEGmods))

write.table(plotDat,file="enrichmentPerModule_EGcell_NLG_v2_TADA_FDR_updated.txt",sep='\t',quote=F,row.names=F)


#myColors=c("darkgray","turquoise3","#FF6666","springgreen3")
png("brainSpanModules_enrichment_EGcell_NEG_v2_TADA_1_largeFont_updated.png",height=7.5,width=10.5,units='in',res=500)
ggplot(plotDat,aes(x= Module, y = FDR, fill=Category, order=Category)) + 
  geom_bar(data=plotDat[which(plotDat$Category=="1EG enrichment" | plotDat$Category=="3NEG enrichment"),] , stat = "identity") + 
  geom_bar(data=plotDat[which(plotDat$Category=="2Potential ASD gene enrichment"),], stat = "identity") +
  scale_fill_manual(name="",values=c("#F8766D","#00BA38","#00BFC4","black")) +
  geom_abline(intercept=-log10(0.1),slope=0, linetype="dashed",color="red") +
  geom_abline(intercept=log10(0.1),slope=0, linetype="dashed",color="red") +
  geom_abline(intercept=0,slope=0, linetype="solid",color="gray50") +
  theme_bw() + xlab("Module ID") + ylab("-log10(q)")+ ggtitle("")+
  theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=18),axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(breaks=rev(c(40,35,30,25,20,15,10,5,1,0,-1,-5)),labels=rev(c("40","35","30","25","20","15","10","5","1","0","1","5")))
  #+ coord_flip()
dev.off()

