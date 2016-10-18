#run in R 3.1.2
library(network)
library(sna)
library(ggplot2)
library(GGally)

setwd("...")
datNode=read.table("Coexp_3modules_EG_nodes_0.8.txt.EGannot",header=T,sep='\t',na.strings=c(""),quote="")
datNode_EG=datNode[which(datNode$Essentialily=="EG"),]
datNode_NEG=datNode[which(datNode$Essentialily=="NEG"),]

degreeDic=datNode$NumberOfUndirectedEdges
names(degreeDic)=datNode$name
geneSetDic=datNode$Essentialily
names(geneSetDic)=datNode$name
geneNameDic=datNode$gene.name
names(geneNameDic)=datNode$name
sfariDic=datNode$SFARI_score
names(sfariDic)=datNode$name
modDic=datNode$Module
names(modDic)=datNode$name
TADAdic=as.numeric(levels(datNode$TADA_FDR))[datNode$TADA_FDR]
names(TADAdic)=datNode$name

wilcox.test(datNode_EG$Degree,datNode_NEG$Degree,alternative="greater") #2.289e-09


datEdge=read.table("Coexp_3modules_EG_edges_0.8.txt",header=T,sep='\t',na.strings=c(""))
edgeList=data.frame(do.call('rbind', strsplit(as.character(datEdge$name),'|',fixed=TRUE)))
edgeList$weights= ifelse(edgeList$X3=="Physical interactions", 1,2)
net=network(edgeList[,c(1,2,4)],directed = FALSE,matrix.type="edgelist",ignore.eval = FALSE, names.eval = "weights")
net %v% "GeneSet" = "Unknown"
net %v% "Degree" = NA
net %v% "Gene" = NA
net %v% "SFARIgene" = NA
net %v% "SFARIscore" = NA
net %v% "Color" = NA
net %v% "Module" = NA
net %v% "TADA" = NA
for (i in 1:length(net$val)){
	geneID=net$val[[i]]$vertex.names
	net$val[[i]]$GeneSet=geneSetDic[geneID]
	net$val[[i]]$Degree=degreeDic[geneID]
	net$val[[i]]$Gene=geneNameDic[geneID]
	net$val[[i]]$Module=modDic[geneID]
	#if (sfariDic[geneID] %in% c("1","2","3","1S","2S","3S","S")){
	if (sfariDic[geneID] %in% c("1","2","1S","2S","S")){
		net$val[[i]]$SFARIgene=paste(geneNameDic[geneID],'(',sfariDic[geneID],')',sep='')
	}
	else {net$val[[i]]$SFARIgene=""}
	if (sfariDic[geneID] %in% c("S","1","1S","2","2S","3","3S","4","4S")){
		net$val[[i]]$SFARIscore=toString(sfariDic[geneID])
	}
	else {net$val[[i]]$SFARIscore=""}
	
	if (geneSetDic[geneID]=="EG"){
		net$val[[i]]$Color="#F8766D"}
	else if (geneSetDic[geneID]=="NEG"){
		net$val[[i]]$Color="#00BFC4"}
	else {net$val[[i]]$Color="grey60"}
	
	if (is.na(TADAdic[geneID])){
		net$val[[i]]$TADA=""}
	else if (as.numeric(TADAdic[geneID])<=0.1){
		net$val[[i]]$TADA="1"}
	else if (as.numeric(TADAdic[geneID])<=0.5){
		net$val[[i]]$TADA="2"}
	else {net$val[[i]]$TADA=""}
}

set.edge.attribute(net,"color", ifelse(net %e% "weights"==1,"lightpink","lightcyan3")) #Physicial interaction and co-expression

#Network plot
png("geneNetwork_coexp_EGcell_NEG_v2_4modules.png",height=6,width=8,units='in',res=800)
ggnet2(net, mode = "fruchtermanreingold", size='degree', max_size=3.5, size.cut=FALSE, color="GeneSet", color.palette=c("EG" = "#F8766D", "NEG" = "#00BFC4", "Unknown"="grey60"), label="Module", label.color="black", label.size=1.5, alpha="GeneSet",alpha.palette = c("EG" = 0.9, "NEG" = 0.8, "Unknown"=0.6), edge.color = "color",edge.size=0.1)
dev.off()

#Network plot:EG only
png("geneNetwork_coexp_0.8_EG_3modules_TADAlabel_test.png",height=6,width=7,units='in',res=800)
ggnet2(net, mode = "fruchtermanreingold", size=2, size.cut=FALSE, color="Module", color.palette=c("M01" = "#E58700", "M02" = "#B983FF", "M16" = "#00BA38"), label="TADA", label.color="black", label.size=1, edge.color = "color",edge.size=0.1,alpha=0.8)
dev.off()

plotDat0=data.frame(datNode_EG$Degree,rep('EG',length(datNode_EG$Degree)))
colnames(plotDat0) <- c("V1","V2")
plotDat1=data.frame(datNode_NEG$Degree,rep('NEG',length(datNode_NEG$Degree)))
colnames(plotDat1) <- c("V1","V2")
plotDat=rbind(plotDat0,plotDat1)
#Degree histogram
png("distDegree_EGcell_NEG_v2_4modules_TADA0.5.png",height=6,width=7,units='in',res=300)
ggplot(plotDat,aes(x=V1,group=V2,fill=as.factor(V2))) + 
  geom_histogram(aes(y = ..density..),position="identity", alpha=0.5,right=TRUE, binwidth=1) + 
  scale_fill_discrete(name="Histogram") +
  geom_line(aes(x=V1,group=V2,color=V2), stat = 'density',position = "identity",adjust=0.5) +
  theme_bw() + xlab("Node degree") + ylab("Density")+ ggtitle("") +
  theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=18))
dev.off()
#Degree boxplot
png("boxplot_nodeDegree_EGcell_NEG_v2.png",height=6,width=4,units='in',res=300)
ggplot(plotDat, aes(x=V2,y=V1,fill=as.factor(V2))) +
	geom_boxplot(notch = TRUE) +
	theme_bw() + xlab("") + ylab("Gene connectivity")+ ggtitle("") +
	theme(plot.title=element_text(size=20),legend.text=element_text(size=14),axis.text=element_text(size=14),axis.title=element_text(size=18))

dev.off()
