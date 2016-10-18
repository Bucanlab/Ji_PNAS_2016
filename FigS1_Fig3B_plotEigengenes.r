library(ggplot2)
library(reshape2)
library(WGCNA);
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

##Load Expression data and gene set

geneSetName="allGenes"
workingDir = paste("...",geneSetName,sep="")
geneSetF = paste(".",geneSetName,".txt",sep="")
dir.create(workingDir,showWarnings = F)
setwd(workingDir)
geneSet = read.table(geneSetF,header=T)
load(paste("...","brainSpan_RNAseq_datExpr3_proteinCoding_fullStages.RData",sep=""))
currentGenes=colnames(datExpr3)
sum(currentGenes %in% geneSet$Gene_symbol)
datExpr4=datExpr3[,currentGenes %in% geneSet$Gene_symbol]
dim(datExpr4)

geneSetName="EG_all_cell_v2"
geneSetF = paste("...",geneSetName,".txt",sep="")
geneSet = read.table(geneSetF,header=T)
datExpr_EG=datExpr3[,currentGenes %in% geneSet$Gene_symbol]
geneSetName="NEG_all_v2"
geneSetF = paste("...",geneSetName,".txt",sep="")
geneSet = read.table(geneSetF,header=T)
datExpr_NEG=datExpr3[,currentGenes %in% geneSet$Gene_symbol]

#write.table(datExpr4,file=paste("detExpr4_",geneSetName,"_fullStages.txt",sep=""),quote=F,row.names=F)
geneSetName="allGenes"


##One-step network construction and module detection
#Caution: determine power from the previous step!
net = blockwiseModules(datExpr4, power = 6, TOMType = "signed", networkType = "signed", maxBlockSize = 40000, saveTOMs = TRUE, saveTOMFileBase = "brainSpan_RNAseq_signed_6_proteinCoding_fullStages", verbose = 3, deepSplit = 2)
sort(table(net$colors),decreasing=T)

moduleLabels = net$colors
# mergedColors = labels2colors(net$colors) # Convert labels to colors for plotting
# Plot the dendrogram and the module colors underneath
png(file=paste("dentrogram_",geneSetName,"_fullStages.png",sep=""),height=9,width=12,units='in',res=500)
p=plotDendroAndColors(net$dendrograms[[1]], moduleLabels[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

MEs = net$MEs;
rownames(MEs)=rownames(datExpr4)
geneTree = net$dendrograms[[1]];
save(net,MEs, moduleLabels, geneTree, file = paste("brainSpan_RNAseq_networkConstruction_auto_",geneSetName,"_fullStages.RData",sep=""))
load(paste("brainSpan_RNAseq_networkConstruction_auto_",geneSetName,"_fullStages.RData",sep=""))

moduleStatL=sort(table(moduleLabels),decreasing=T)
moduleNameL=names(moduleStatL)
moduleLabels=factor(moduleLabels)
#Sort and rename modules
for (moduleIdx in 1:length(moduleNameL)){
	eigenGene=moduleNameL[moduleIdx]
	oldName=paste(eigenGene,sep="")
	if (moduleIdx<10){
		newName=paste("M0",moduleIdx,"_",eigenGene,"_",moduleStatL[moduleIdx],sep="")}
	else {newName=paste("M",moduleIdx,"_",eigenGene,"_",moduleStatL[moduleIdx],sep="") }
	levels(moduleLabels)[which(levels(moduleLabels)==oldName)] <- newName
}


#Output gene module membership
membershipDf=data.frame(colnames(datExpr4),moduleLabels)
colnames(membershipDf)=c("Gene","Module")
write.table(membershipDf,file=paste("moduleMembership_",geneSetName,"_fullStages.txt",sep=""),quote=F,row.names=F)

#load(file = paste("brainSpan_RNAseq_networkConstruction_auto_",geneSetName,"_fullStages.RData",sep=""))
#load(file = "brainSpan_RNAseq_signed_6_proteinCoding_fullStages-block.1.RData")

##BrainSpan Figures
#read Brainspan column info
tissueDat=read.table(".../columns_metadata.csv",header=T,sep=',')
MEs$tissueID=rownames(MEs)
MEs$age=tissueDat[MEs$tissueID,]$age
MEs$tissue=tissueDat[MEs$tissueID,]$structure_name
MEs_long=melt(MEs)
colnames(MEs_long)[4]="eigenGene"
colnames(MEs_long)[5]="PC1"
tissueOrder=c("8 pcw","9 pcw","12 pcw","13 pcw","16 pcw","17 pcw","19 pcw","21 pcw","24 pcw","25 pcw","26 pcw","35 pcw","37 pcw","4 mos","10 mos","1 yrs","2 yrs","3 yrs","4 yrs","8 yrs","11 yrs","13 yrs","15 yrs","18 yrs","19 yrs","21 yrs","23 yrs","30 yrs","36 yrs","37 yrs","40 yrs")
MEs_long$age=factor(MEs_long$age,tissueOrder)

for (moduleIdx in 1:length(moduleNameL)){
	eigenGene=moduleNameL[moduleIdx]
	oldName=paste("ME",eigenGene,sep="")
	#if (moduleIdx<10){
		#newName=paste("M0",moduleIdx,"_",eigenGene,"_",moduleStatL[moduleIdx],sep="")}
	#else {newName=paste("M",moduleIdx,"_",eigenGene,"_",moduleStatL[moduleIdx],sep="") }
	if (moduleIdx<10){
		newName=paste("M0",moduleIdx," (n=",moduleStatL[moduleIdx],")",sep="") }
	else {newName=paste("M",moduleIdx," (n=",moduleStatL[moduleIdx],")",sep="") }
	levels(MEs_long$eigenGene)[which(levels(MEs_long$eigenGene)==oldName)] <- newName
}

MEs_long$eigenGene=factor(MEs_long$eigenGene,levels=sort(levels(MEs_long$eigenGene)))

#Facet: all figures in one 
png(paste("MEs_",geneSetName,"_fullstages","_loess_largeFont.png",sep=""),height=10,width=16,units='in',res=500)
p = ggplot(MEs_long,aes(x=age,y=PC1,colour=factor(tissue))) +
	geom_point(size=1.0) +
	scale_color_discrete("") +
	geom_smooth(method="loess",aes(group=tissue),se=FALSE) +  #Polynomial: method="lm",formula=y~poly(x,3)
	theme_bw() + 
	theme(axis.text.x = element_blank(),axis.title.x=element_blank()) +
	labs(aesthetic="") +
	geom_vline(xintercept=13.5,linetype = "dashed",colour="black") +
	ylim(-0.1,0.3) +
	facet_wrap(~ eigenGene,scales='free') +
	theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0), legend.position="right", legend.title = element_text(size=12),legend.text = element_text(size=12))
print(p)
dev.off()

#One figure per module
for (moduleIdx in 1:length(moduleNameL)){
	eigenGene=moduleNameL[moduleIdx]
	MEs_long_this=MEs_long[which(grepl(paste('_',eigenGene,'_',sep=''),MEs_long$eigenGene)),]
	png(paste("MEs_M",moduleIdx,'_',moduleStatL[moduleIdx],'_',eigenGene,"_loess.png",sep=""),height=9,width=15,units='in',res=500)
	p = ggplot(MEs_long_this,aes(x=age,y=PC1,colour=factor(tissue))) +
		geom_point() +
		scale_color_discrete("Brain tissue") +
		geom_smooth(method="loess",aes(group=tissue),se=FALSE) +  #Polynomial: method="lm",formula=y~poly(x,3)
		labs(aesthetic='Brain tissue') +
		geom_vline(xintercept=13.5,linetype = "dashed",colour="red") +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
		theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=18))
	print(p)
	dev.off()
}

geneSetName="allGenes"
# 3 modules in one figure:Pooled tissues, PC1 from allGenes
MEs_long_subset=MEs_long[which(MEs_long$eigenGene %in% c("M01_turquoise_1601","M02_blue_1150","M16_midnightblue_347")),]
png(paste("MEs_",geneSetName,"_fullstages","_loess_pooled_3modules_PC1_allGenes.png",sep=""),height=6,width=9,units='in',res=500)
p = ggplot(MEs_long_subset,aes(x=age,y=PC1,colour=factor(eigenGene))) +
	#geom_point() +
	scale_color_discrete("Module") +
	geom_smooth(method="loess",aes(group=eigenGene),se=FALSE) +  #Polynomial: method="lm",formula=y~poly(x,3)
	labs(aesthetic='Module') +
	geom_vline(xintercept=13.5,linetype = "dashed",colour="black") +
	xlab("") + ylab("1st PC of Expression Profiles") + ylim(-0.05,0.10) +
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
	theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=18))
print(p)
dev.off()

##Recalculate PC1 for EGs and NEGs seperately
geneMembership=read.table("moduleMembership_allGenes_fullStages.txt.EGannot",header=T,sep="\t")
#Recalculate PC1 for EGs
memVector_EG=geneMembership$Module[match(colnames(datExpr_EG),geneMembership$Gene)]
MEs_EG_raw=moduleEigengenes(datExpr_EG, memVector_EG)
MEs_EG=MEs_EG_raw$eigengenes
MEs_EG$tissueID=rownames(MEs_EG)
MEs_EG$age=tissueDat[MEs_EG$tissueID,]$age
MEs_EG$tissue=tissueDat[MEs_EG$tissueID,]$structure_name
MEs_EG_long=melt(MEs_EG)
colnames(MEs_EG_long)[4]="eigenGene"
colnames(MEs_EG_long)[5]="PC1"
tissueOrder=c("8 pcw","9 pcw","12 pcw","13 pcw","16 pcw","17 pcw","19 pcw","21 pcw","24 pcw","25 pcw","26 pcw","35 pcw","37 pcw","4 mos","10 mos","1 yrs","2 yrs","3 yrs","4 yrs","8 yrs","11 yrs","13 yrs","15 yrs","18 yrs","19 yrs","21 yrs","23 yrs","30 yrs","36 yrs","37 yrs","40 yrs")
MEs_EG_long$age=factor(MEs_EG_long$age,tissueOrder)
MEs_EG_long$eigenGene=factor(MEs_EG_long$eigenGene,levels=sort(levels(MEs_EG_long$eigenGene)))
# 4 modules in one figure:Pooled tissues, PC1 from EGs
MEs_EG_long_subset=MEs_EG_long[which(MEs_EG_long$eigenGene %in% c("MEM01_turquoise_1601","MEM02_blue_1150","MEM16_midnightblue_347")),]
#Three modules in one plot
png(paste("MEs_",geneSetName,"_fullstages","_loess_pooled_3modules_PC1_EG.png",sep=""),height=6,width=9,units='in',res=500)
p = ggplot(MEs_EG_long_subset,aes(x=age,y=PC1,colour=factor(eigenGene))) +
	#geom_point() +
	#scale_color_discrete("Module") +
	scale_colour_manual(values=c("#E58700", "#B983FF", "#00BA38")) +
	geom_smooth(method="loess",aes(group=eigenGene),se=FALSE) +  #Polynomial: method="lm",formula=y~poly(x,3)
	labs(aesthetic='Module') +
	geom_vline(xintercept=13.5,linetype = "dashed",colour="black") +
	xlab("") + ylab("1st PC of Expression Profiles") + ylim(-0.05,0.10) +
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
	theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=18))
print(p)
dev.off()

#Facet: all figures in one: PC1 for EGs
png(paste("MEs_",geneSetName,"_fullstages","_loess_PC1_EGs.png",sep=""),height=12,width=18,units='in',res=500)
p = ggplot(MEs_EG_long,aes(x=age,y=PC1,colour="black")) +
	geom_smooth(method="loess",aes(group=eigenGene),se=FALSE) +  #Polynomial: method="lm",formula=y~poly(x,3)
	theme_bw() + 
	theme(axis.text.x = element_blank(),axis.title.x=element_blank()) +
	geom_vline(xintercept=13.5,linetype = "dashed",colour="red") +
	ylim(-0.1,0.3) +
	facet_wrap(~ eigenGene,scales='free')
print(p)
dev.off()

#Recalculate PC1 for NEGs
memVector_NEG=geneMembership$Module[match(colnames(datExpr_NEG),geneMembership$Gene)]
MEs_NEG_raw=moduleEigengenes(datExpr_NEG, memVector_NEG)
MEs_NEG=MEs_NEG_raw$eigengenes
MEs_NEG$tissueID=rownames(MEs_NEG)
MEs_NEG$age=tissueDat[MEs_NEG$tissueID,]$age
MEs_NEG$tissue=tissueDat[MEs_NEG$tissueID,]$structure_name
MEs_NEG_long=melt(MEs_NEG)
colnames(MEs_NEG_long)[4]="eigenGene"
colnames(MEs_NEG_long)[5]="PC1"
tissueOrder=c("8 pcw","9 pcw","12 pcw","13 pcw","16 pcw","17 pcw","19 pcw","21 pcw","24 pcw","25 pcw","26 pcw","35 pcw","37 pcw","4 mos","10 mos","1 yrs","2 yrs","3 yrs","4 yrs","8 yrs","11 yrs","13 yrs","15 yrs","18 yrs","19 yrs","21 yrs","23 yrs","30 yrs","36 yrs","37 yrs","40 yrs")
MEs_NEG_long$age=factor(MEs_NEG_long$age,tissueOrder)
MEs_NEG_long$eigenGene=factor(MEs_NEG_long$eigenGene,levels=sort(levels(MEs_NEG_long$eigenGene)))
# 4 modules in one figure:Pooled tissues, PC1 from NEGs
MEs_NEG_long_subset=MEs_NEG_long[which(MEs_NEG_long$eigenGene %in% c("MEM01_turquoise_1601","MEM02_blue_1150","MEM16_midnightblue_347","MEM17_lightcyan_339")),]
png(paste("MEs_",geneSetName,"_fullstages","_loess_pooled_4modules_PC1_NEG.png",sep=""),height=6,width=9,units='in',res=500)
p = ggplot(MEs_NEG_long_subset,aes(x=age,y=PC1,colour=factor(eigenGene))) +
	#geom_point() +
	scale_color_discrete("Module") +
	geom_smooth(method="loess",aes(group=eigenGene),se=FALSE) +  #Polynomial: method="lm",formula=y~poly(x,3)
	labs(aesthetic='Module') +
	geom_vline(xintercept=13.5,linetype = "dashed",colour="black") +
	xlab("") + ylab("1st PC of Expression Profiles") + ylim(-0.05,0.10) +
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
	theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=18))	
print(p)
dev.off()


##GTEx Figures
#prepare dataframe for ggplot
MEs$tissue=rownames(MEs)
MEs_long=melt(MEs)
colnames(MEs_long)[2]="eigenGene"
colnames(MEs_long)[3]="PC1"

for (moduleIdx in 1:length(moduleNameL)){
	eigenGene=moduleNameL[moduleIdx]
	oldName=paste("ME",eigenGene,sep="")
	if (moduleIdx<10){
		newName=paste("M0",moduleIdx,"_",eigenGene,"_",moduleStatL[moduleIdx],sep="")}
	else {newName=paste("M",moduleIdx,"_",eigenGene,"_",moduleStatL[moduleIdx],sep="") }
	levels(MEs_long$eigenGene)[which(levels(MEs_long$eigenGene)==oldName)] <- newName
}

MEs_long$eigenGene=factor(MEs_long$eigenGene,levels=sort(levels(MEs_long$eigenGene)))

#Facet: all figures in one 
png(paste("MEs_",geneSetName,"_GTEx.png",sep=""),height=12,width=18,units='in',res=500)
p = ggplot(MEs_long,aes(x=tissue,y=PC1,fill=factor(tissue))) +
	geom_bar(stat="identity") +
	scale_color_discrete("Tissue") +
	theme_bw() + 
	theme(axis.text.x = element_blank(),axis.title.x=element_blank()) +
	labs(aesthetic='Tissue') +
	#ylim(-0.1,0.3) +
	facet_wrap(~ eigenGene,scales='free')
print(p)
dev.off()

#One figure per module
for (moduleIdx in 1:length(moduleNameL)){
	eigenGene=moduleNameL[moduleIdx]
	MEs_long_this=MEs_long[which(grepl(eigenGene,MEs_long$eigenGene)),]
	png(paste("MEs_M",moduleIdx,'_',moduleStatL[moduleIdx],'_',eigenGene,".png",sep=""),height=9,width=15,units='in',res=300)
	p = ggplot(MEs_long_this,aes(x=age,y=PC1,colour=factor(tissue))) +
		geom_point() +
		scale_color_discrete("Brain tissue") +
		geom_smooth(method="loess",aes(group=tissue),se=FALSE) +  #Polynomial: method="lm",formula=y~poly(x,3)
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
		labs(aesthetic='Brain tissue') +
		geom_vline(xintercept=13.5,linetype = "dashed",colour="red")
	print(p)
	dev.off()
}



#Output hard-thresholding co-expression network
library(reshape2)
adjMat=adjacency(datExpr4,selectCols = NULL, type = "signed", power = 1, corFnc = "cor")
threshold=0.99
adjMat_th=signumAdjacencyFunction(adjMat,threshold)
colnames(adjMat_th) =  names(datExpr4)
rownames(adjMat_th) =  names(datExpr4)
adjMat_th[adjMat_th==0]=NA
adjMat_th_out=melt(adjMat_th, na.rm=T)
write.table(adjMat_th_out,file=paste("brainSpan_RNAseq_adjacency_hardTH_",threshold,"_proteinCoding.txt",sep=""),quote=F,row.names=F)

