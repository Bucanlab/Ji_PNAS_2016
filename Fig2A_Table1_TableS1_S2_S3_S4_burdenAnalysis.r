library(ggplot2)
library(lsr)

bootstrapCohensD_paired <- function(vector1,vector2,n){
	cohensDVec=c()
	vecLen=length(vector1)
	
	for (i in 1:n){
		#print(i)
		thisSample=sample(c(1:vecLen),replace=T)
		thisCohensD=cohensD(vector1[thisSample],vector2[thisSample],method = "paired")
		mean1=mean(vector1[thisSample])
		mean2=mean(vector2[thisSample])
		if (mean1<mean2){
			thisCohensD=-thisCohensD
		}
		cohensDVec=append(cohensDVec,thisCohensD)
	}
	cohensDVec=sort(cohensDVec)
	return(cohensDVec)
}

setwd("...")
dat=read.table("masterPheno_indivBurden_final.txt",header=T,sep='\t')
dat_quad=dat[which(dat$Family.Type=="quad"),]

###Figure 2A, Table S1.
## Burden difference between pairs of siblings: EG
dim(dat_quad) #1781
mean(dat_quad$Burden_proband_EG_dnLoF) #0.06400898
mean(dat_quad$Burden_sibling_EG_dnLoF) #0.0286356
sd(dat_quad$Burden_proband_EG_dnLoF) #0.2538501
sd(dat_quad$Burden_sibling_EG_dnLoF) #0.1701612
cohensD(dat_quad$Burden_proband_EG_dnLoF,dat_quad$Burden_sibling_EG_dnLoF,method = "paired") #0.1169762
cohensDVec=bootstrapCohensD_paired(dat_quad$Burden_proband_EG_dnLoF,dat_quad$Burden_sibling_EG_dnLoF,10000)
cohensDVec[251] #0.07154123
cohensDVec[9750] #0.1595998
wilcox.test(dat_quad$Burden_proband_EG_dnLoF, dat_quad$Burden_sibling_EG_dnLoF, paired=TRUE, alternative="greater") #p-value = 4.751e-07

dim(dat_quad) #1781
mean(dat_quad$Burden_proband_EG_dnNSD) #0.206064
mean(dat_quad$Burden_sibling_EG_dnNSD) #0.1588995
cohensD(dat_quad$Burden_proband_EG_dnNSD,dat_quad$Burden_sibling_EG_dnNSD,method="paired") #0.07940217
cohensDVec=bootstrapCohensD_paired(dat_quad$Burden_proband_EG_dnNSD,dat_quad$Burden_sibling_EG_dnNSD,10000)
cohensDVec[251] #0.03238893
cohensDVec[9750] #0.1273962
wilcox.test(dat_quad$Burden_proband_EG_dnNSD, dat_quad$Burden_sibling_EG_dnNSD, paired=TRUE, alternative="greater") #p-value = 0.0003405

dim(dat_quad) #1781
mean(dat_quad$Burden_proband_EG_inhRD) #10.74284
mean(dat_quad$Burden_sibling_EG_inhRD) #10.60415
cohensD(dat_quad$Burden_proband_EG_inhRD,dat_quad$Burden_sibling_EG_inhRD,method="paired") #0.04199952
cohensDVec=bootstrapCohensD_paired(dat_quad$Burden_proband_EG_inhRD,dat_quad$Burden_sibling_EG_inhRD,10000)
cohensDVec[251] #-0.00410884
cohensDVec[9750] #0.08866722
wilcox.test(dat_quad$Burden_proband_EG_inhRD, dat_quad$Burden_sibling_EG_inhRD, paired=TRUE, alternative="greater") #p-value = 0.01688


#Burden difference between pairs of siblings: NEG
dim(dat_quad) #1781
mean(dat_quad$Burden_proband_NEG_dnLoF) #0.03874228
mean(dat_quad$Burden_sibling_NEG_dnLoF) #0.03088153
cohensD(dat_quad$Burden_proband_NEG_dnLoF,dat_quad$Burden_sibling_NEG_dnLoF,method = "paired") #0.03003932
cohensDVec=bootstrapCohensD_paired(dat_quad$Burden_proband_NEG_dnLoF,dat_quad$Burden_sibling_NEG_dnLoF,10000)
cohensDVec[251] #-0.01574115
cohensDVec[9750] #0.07744043
wilcox.test(dat_quad$Burden_proband_NEG_dnLoF, dat_quad$Burden_sibling_NEG_dnLoF, paired=TRUE, alternative="greater") #p-value = 0.1028

dim(dat_quad) #1781
mean(dat_quad$Burden_proband_NEG_dnNSD) #0.1611454
mean(dat_quad$Burden_sibling_NEG_dnNSD) #0.1403706
cohensD(dat_quad$Burden_proband_NEG_dnNSD,dat_quad$Burden_sibling_NEG_dnNSD,method="paired") #0.03743397
cohensDVec=bootstrapCohensD_paired(dat_quad$Burden_proband_NEG_dnNSD,dat_quad$Burden_sibling_NEG_dnNSD,10000)
cohensDVec[251] #-0.01002884
cohensDVec[9750] #0.08266127
wilcox.test(dat_quad$Burden_proband_NEG_dnNSD, dat_quad$Burden_sibling_NEG_dnNSD, paired=TRUE, alternative="greater") #p-value = 0.0691

dim(dat_quad) #1781
mean(dat_quad$Burden_proband_NEG_inhRD) #12.78158
mean(dat_quad$Burden_sibling_NEG_inhRD) #12.83549
cohensD(dat_quad$Burden_proband_NEG_inhRD,dat_quad$Burden_sibling_NEG_inhRD,method="paired") #-0.01503358
cohensDVec=bootstrapCohensD_paired(dat_quad$Burden_proband_NEG_inhRD,dat_quad$Burden_sibling_NEG_inhRD,10000)
cohensDVec[251] #-0.06180866
cohensDVec[9750] #0.03079507
wilcox.test(dat_quad$Burden_proband_NEG_inhRD, dat_quad$Burden_sibling_NEG_inhRD, paired=TRUE, alternative="greater") #p-value = 0.7456

###Table S2
## Burden difference in probands between genders: dnLoF in EG
dat_M_pro=dat[which(dat$probandGender=="M"),]
dat_F_pro=dat[which(dat$probandGender=="F"),]
mean(dat_M_pro$Burden_proband_EG_dnLoF) #0.0597161
mean(dat_F_pro$Burden_proband_EG_dnLoF) #0.08615385
cohensD(dat_F_pro$Burden_proband_EG_dnLoF, dat_M_pro$Burden_proband_EG_dnLoF) #0.1042026
wilcox.test(dat_F_pro$Burden_proband_EG_dnLoF, dat_M_pro$Burden_proband_EG_dnLoF, paired=FALSE, alternative="greater") #p-value = 0.03548

## Burden difference in probands between genders: dnLoF in NEG
mean(dat_M_pro$Burden_proband_NEG_dnLoF) #0.03573177
mean(dat_F_pro$Burden_proband_NEG_dnLoF) #0.04615385
cohensD(dat_F_pro$Burden_proband_NEG_dnLoF, dat_M_pro$Burden_proband_NEG_dnLoF) #0.05508355
wilcox.test(dat_F_pro$Burden_proband_NEG_dnLoF, dat_M_pro$Burden_proband_NEG_dnLoF, paired=FALSE, alternative="greater") #p-value = 0.1782

## Burden difference in probands between genders: dnNSD in EG
mean(dat_M_pro$Burden_proband_EG_dnNSD) #0.1948116
mean(dat_F_pro$Burden_proband_EG_dnNSD) #0.24
cohensD(dat_F_pro$Burden_proband_EG_dnNSD, dat_M_pro$Burden_proband_EG_dnNSD) #0.1014154
wilcox.test(dat_F_pro$Burden_proband_EG_dnNSD, dat_M_pro$Burden_proband_EG_dnNSD, paired=FALSE, alternative="greater") #p-value = 0.03883

## Burden difference in probands between genders: dnNSD in NEG
mean(dat_M_pro$Burden_proband_NEG_dnNSD) #0.1595693
mean(dat_F_pro$Burden_proband_NEG_dnNSD) #0.2
cohensD(dat_F_pro$Burden_proband_NEG_dnNSD, dat_M_pro$Burden_proband_NEG_dnNSD) #0.09933147
wilcox.test(dat_F_pro$Burden_proband_NEG_dnNSD, dat_M_pro$Burden_proband_NEG_dnNSD, paired=FALSE, alternative="greater") #p-value = 0.07416

## Burden difference in probands between genders: inhRD in EG
mean(dat_M_pro$Burden_proband_EG_inhRD) #10.96329
mean(dat_F_pro$Burden_proband_EG_inhRD) #11.05231
cohensD(dat_F_pro$Burden_proband_EG_inhRD, dat_M_pro$Burden_proband_EG_inhRD) #0.01509099
wilcox.test(dat_F_pro$Burden_proband_EG_inhRD, dat_M_pro$Burden_proband_EG_inhRD, paired=FALSE, alternative="greater") #p-value = 0.4711

## Burden difference in probands between genders: inhRD in NEG
mean(dat_M_pro$Burden_proband_NEG_inhRD) #13.01126
mean(dat_F_pro$Burden_proband_NEG_inhRD) #13.26769
cohensD(dat_F_pro$Burden_proband_NEG_inhRD, dat_M_pro$Burden_proband_NEG_inhRD) #0.03596789
wilcox.test(dat_F_pro$Burden_proband_NEG_inhRD, dat_M_pro$Burden_proband_NEG_inhRD, paired=FALSE, alternative="greater") #p-value = 0.5271

###Table S3
## Burden difference between pairs of siblings: dnLoF in EG
dat_quad_FM=dat_quad[which(dat_quad$probandGender=="F" & dat_quad$siblingGender=="M"),]
dat_quad_MF=dat_quad[which(dat_quad$probandGender=="M" & dat_quad$siblingGender=="F"),]
dat_quad_MM=dat_quad[which(dat_quad$probandGender=="M" & dat_quad$siblingGender=="M"),]
dat_quad_FF=dat_quad[which(dat_quad$probandGender=="F" & dat_quad$siblingGender=="F"),]

dim(dat_quad) #1781
mean(dat_quad$Burden_proband_EG_dnLoF) #0.06400898
mean(dat_quad$Burden_sibling_EG_dnLoF) #0.0286356
cohensD(dat_quad$Burden_proband_EG_dnLoF,dat_quad$Burden_sibling_EG_dnLoF,method = "paired") #0.1169762
wilcox.test(dat_quad$Burden_proband_EG_dnLoF, dat_quad$Burden_sibling_EG_dnLoF, paired=TRUE, alternative="greater") #p-value = 6.541e-07

dim(dat_quad_FM) #101
mean(dat_quad_FM$Burden_proband_EG_dnLoF) #0.08910891
mean(dat_quad_FM$Burden_sibling_EG_dnLoF) #0.00990099
cohensD(dat_quad_FM$Burden_proband_EG_dnLoF,dat_quad_FM$Burden_sibling_EG_dnLoF,method = "paired") #0.2588116
wilcox.test(dat_quad_FM$Burden_proband_EG_dnLoF, dat_quad_FM$Burden_sibling_EG_dnLoF, paired=TRUE, alternative="greater") #p-value = 0.006712

dim(dat_quad_MF) #826
mean(dat_quad_MF$Burden_proband_EG_dnLoF) #0.05932203
mean(dat_quad_MF$Burden_sibling_EG_dnLoF) #0.03268765
cohensD(dat_quad_MF$Burden_proband_EG_dnLoF,dat_quad_MF$Burden_sibling_EG_dnLoF,method = "paired") #0.08928531
wilcox.test(dat_quad_MF$Burden_proband_EG_dnLoF, dat_quad_MF$Burden_sibling_EG_dnLoF, paired=TRUE, alternative="greater") #p-value = 0.005319

dim(dat_quad_MM) #732
mean(dat_quad_MM$Burden_proband_EG_dnLoF) #0.06147541
mean(dat_quad_MM$Burden_sibling_EG_dnLoF) #0.02459016
cohensD(dat_quad_MM$Burden_proband_EG_dnLoF,dat_quad_MM$Burden_sibling_EG_dnLoF,method = "paired") #0.1227512
wilcox.test(dat_quad_MM$Burden_proband_EG_dnLoF, dat_quad_MM$Burden_sibling_EG_dnLoF, paired=TRUE, alternative="greater") #p-value = 0.0005097

dim(dat_quad_FF) #122
mean(dat_quad_FF$Burden_proband_EG_dnLoF) #0.09016393
mean(dat_quad_FF$Burden_sibling_EG_dnLoF) #0.04098361
cohensD(dat_quad_FF$Burden_proband_EG_dnLoF,dat_quad_FF$Burden_sibling_EG_dnLoF,method = "paired") #0.1461322
wilcox.test(dat_quad_FF$Burden_proband_EG_dnLoF, dat_quad_FF$Burden_sibling_EG_dnLoF, paired=TRUE, alternative="greater") #p-value = 0.05998


## Burden difference between pairs of siblings: dnLoF in NEG
dim(dat_quad) #1781
mean(dat_quad$Burden_proband_NEG_dnLoF) #0.03874228
mean(dat_quad$Burden_sibling_NEG_dnLoF) #0.03088153
cohensD(dat_quad$Burden_proband_NEG_dnLoF,dat_quad$Burden_sibling_NEG_dnLoF,method = "paired") #0.03003932
wilcox.test(dat_quad$Burden_proband_NEG_dnLoF, dat_quad$Burden_sibling_NEG_dnLoF, paired=TRUE, alternative="greater") #p-value = 0.09182

dim(dat_quad_FM) #101
mean(dat_quad_FM$Burden_proband_NEG_dnLoF) #0.03960396
mean(dat_quad_FM$Burden_sibling_NEG_dnLoF) #0.02970297
cohensD(dat_quad_FM$Burden_proband_NEG_dnLoF,dat_quad_FM$Burden_sibling_NEG_dnLoF,method = "paired") #0.03744872
wilcox.test(dat_quad_FM$Burden_proband_NEG_dnLoF, dat_quad_FM$Burden_sibling_NEG_dnLoF, paired=TRUE, alternative="greater") #p-value = 0.3884

dim(dat_quad_MF) #826
mean(dat_quad_MF$Burden_proband_NEG_dnLoF) #0.04116223
mean(dat_quad_MF$Burden_sibling_NEG_dnLoF) #0.02663438
cohensD(dat_quad_MF$Burden_proband_NEG_dnLoF,dat_quad_MF$Burden_sibling_NEG_dnLoF,method = "paired") #0.05584846
wilcox.test(dat_quad_MF$Burden_proband_NEG_dnLoF, dat_quad_MF$Burden_sibling_NEG_dnLoF, paired=TRUE, alternative="greater") #p-value = 0.05492

dim(dat_quad_MM) #732
mean(dat_quad_MM$Burden_proband_NEG_dnLoF) #0.03688525
mean(dat_quad_MM$Burden_sibling_NEG_dnLoF) #0.03415301
cohensD(dat_quad_MM$Burden_proband_NEG_dnLoF,dat_quad_MM$Burden_sibling_NEG_dnLoF,method = "paired") #0.02132976
wilcox.test(dat_quad_MM$Burden_proband_NEG_dnLoF, dat_quad_MM$Burden_sibling_NEG_dnLoF, paired=TRUE, alternative="greater") #p-value = 0.2838

dim(dat_quad_FF) #122
mean(dat_quad_FF$Burden_proband_NEG_dnLoF) #0.03278689
mean(dat_quad_FF$Burden_sibling_NEG_dnLoF) #0.05737705
cohensD(dat_quad_FF$Burden_proband_NEG_dnLoF,dat_quad_FF$Burden_sibling_NEG_dnLoF,method = "paired") #0.08183121
wilcox.test(dat_quad_FF$Burden_proband_NEG_dnLoF, dat_quad_FF$Burden_sibling_NEG_dnLoF, paired=TRUE, alternative="greater") #p-value = 0.8302

## Burden difference between pairs of siblings: dnNSD in EG

dim(dat_quad) #1781
mean(dat_quad$Burden_proband_EG_dnNSD) #0.206064
mean(dat_quad$Burden_sibling_EG_dnNSD) #0.1588995
cohensD(dat_quad$Burden_proband_EG_dnNSD,dat_quad$Burden_sibling_EG_dnNSD,method = "paired") #0.07940217
wilcox.test(dat_quad$Burden_proband_EG_dnNSD, dat_quad$Burden_sibling_EG_dnNSD, paired=TRUE, alternative="greater") #p-value = 0.0003405

dim(dat_quad_FM) #101
mean(dat_quad_FM$Burden_proband_EG_dnNSD) #0.2178218
mean(dat_quad_FM$Burden_sibling_EG_dnNSD) #0.1683168
cohensD(dat_quad_FM$Burden_proband_EG_dnNSD,dat_quad_FM$Burden_sibling_EG_dnNSD,method = "paired") #0.07240129
wilcox.test(dat_quad_FM$Burden_proband_EG_dnNSD, dat_quad_FM$Burden_sibling_EG_dnNSD, paired=TRUE, alternative="greater") #p-value = 0.2392

dim(dat_quad_MF) #826
mean(dat_quad_MF$Burden_proband_EG_dnNSD) #0.2094431
mean(dat_quad_MF$Burden_sibling_EG_dnNSD) #0.1755448
cohensD(dat_quad_MF$Burden_proband_EG_dnNSD,dat_quad_MF$Burden_sibling_EG_dnNSD,method = "paired") # 0.05520638
wilcox.test(dat_quad_MF$Burden_proband_EG_dnNSD, dat_quad_MF$Burden_sibling_EG_dnNSD, paired=TRUE, alternative="greater") #p-value = 0.04544

dim(dat_quad_MM) #732
mean(dat_quad_MM$Burden_proband_EG_dnNSD) #0.1885246
mean(dat_quad_MM$Burden_sibling_EG_dnNSD) #0.1270492
cohensD(dat_quad_MM$Burden_proband_EG_dnNSD,dat_quad_MM$Burden_sibling_EG_dnNSD,method = "paired") #0.1135576
wilcox.test(dat_quad_MM$Burden_proband_EG_dnNSD, dat_quad_MM$Burden_sibling_EG_dnNSD, paired=TRUE, alternative="greater") #p-value = 0.001304

dim(dat_quad_FF) #122
mean(dat_quad_FF$Burden_proband_EG_dnNSD) #0.2786885
mean(dat_quad_FF$Burden_sibling_EG_dnNSD) #0.2295082
cohensD(dat_quad_FF$Burden_proband_EG_dnNSD,dat_quad_FF$Burden_sibling_EG_dnNSD,method = "paired") #0.0724832
wilcox.test(dat_quad_FF$Burden_proband_EG_dnNSD, dat_quad_FF$Burden_sibling_EG_dnNSD, paired=TRUE, alternative="greater") #p-value = 0.2157

## Burden difference between pairs of siblings: dnNSD in NEG
dim(dat_quad) #1781
mean(dat_quad$Burden_proband_NEG_dnNSD) #0.1611454
mean(dat_quad$Burden_sibling_NEG_dnNSD) #0.1403706
cohensD(dat_quad$Burden_proband_NEG_dnNSD,dat_quad$Burden_sibling_NEG_dnNSD,method = "paired") #0.03743397
wilcox.test(dat_quad$Burden_proband_NEG_dnNSD, dat_quad$Burden_sibling_NEG_dnNSD, paired=TRUE, alternative="greater") #p-value = 0.0691

dim(dat_quad_FM) #101
mean(dat_quad_FM$Burden_proband_NEG_dnNSD) #0.1881188
mean(dat_quad_FM$Burden_sibling_NEG_dnNSD) #0.1980198
cohensD(dat_quad_FM$Burden_proband_NEG_dnNSD,dat_quad_FM$Burden_sibling_NEG_dnNSD,method = "paired") #0.01546462
wilcox.test(dat_quad_FM$Burden_proband_NEG_dnNSD, dat_quad_FM$Burden_sibling_NEG_dnNSD, paired=TRUE, alternative="greater") #p-value = 0.3884

dim(dat_quad_MF) #826
mean(dat_quad_MF$Burden_proband_NEG_dnNSD) #0.1464891
mean(dat_quad_MF$Burden_sibling_NEG_dnNSD) #0.1476998
cohensD(dat_quad_MF$Burden_proband_NEG_dnNSD,dat_quad_MF$Burden_sibling_NEG_dnNSD,method = "paired") #0.002249309
wilcox.test(dat_quad_MF$Burden_proband_NEG_dnNSD, dat_quad_MF$Burden_sibling_NEG_dnNSD, paired=TRUE, alternative="greater") #p-value = 0.5515

dim(dat_quad_MM) #732
mean(dat_quad_MM$Burden_proband_NEG_dnNSD) #0.1666667
mean(dat_quad_MM$Burden_sibling_NEG_dnNSD) #0.1174863
cohensD(dat_quad_MM$Burden_proband_NEG_dnNSD,dat_quad_MM$Burden_sibling_NEG_dnNSD,method = "paired") #0.09042586
wilcox.test(dat_quad_MM$Burden_proband_NEG_dnNSD, dat_quad_MM$Burden_sibling_NEG_dnNSD, paired=TRUE, alternative="greater") #p-value = 0.008002

dim(dat_quad_FF) #122
mean(dat_quad_FF$Burden_proband_NEG_dnNSD) #0.204918
mean(dat_quad_FF$Burden_sibling_NEG_dnNSD) #0.1803279
cohensD(dat_quad_FF$Burden_proband_NEG_dnNSD,dat_quad_FF$Burden_sibling_NEG_dnNSD,method = "paired") #0.03790385
wilcox.test(dat_quad_FF$Burden_proband_NEG_dnNSD, dat_quad_FF$Burden_sibling_NEG_dnNSD, paired=TRUE, alternative="greater") #p-value = 0.3817

## Burden difference between pairs of siblings: inhRD in EG

dim(dat_quad) #1781
mean(dat_quad$Burden_proband_EG_inhRD) 
mean(dat_quad$Burden_sibling_EG_inhRD) 
cohensD(dat_quad$Burden_proband_EG_inhRD,dat_quad$Burden_sibling_EG_inhRD,method = "paired") 
wilcox.test(dat_quad$Burden_proband_EG_inhRD, dat_quad$Burden_sibling_EG_inhRD, paired=TRUE, alternative="greater") 

dim(dat_quad_FM) #101
mean(dat_quad_FM$Burden_proband_EG_inhRD) 
mean(dat_quad_FM$Burden_sibling_EG_inhRD) 
cohensD(dat_quad_FM$Burden_proband_EG_inhRD,dat_quad_FM$Burden_sibling_EG_inhRD,method = "paired") 
wilcox.test(dat_quad_FM$Burden_proband_EG_inhRD, dat_quad_FM$Burden_sibling_EG_inhRD, paired=TRUE, alternative="greater") 

dim(dat_quad_MF) 
mean(dat_quad_MF$Burden_proband_EG_inhRD) 
mean(dat_quad_MF$Burden_sibling_EG_inhRD) 
cohensD(dat_quad_MF$Burden_proband_EG_inhRD,dat_quad_MF$Burden_sibling_EG_inhRD,method = "paired") 
wilcox.test(dat_quad_MF$Burden_proband_EG_inhRD, dat_quad_MF$Burden_sibling_EG_inhRD, paired=TRUE, alternative="greater") 

dim(dat_quad_MM) 
mean(dat_quad_MM$Burden_proband_EG_inhRD) 
mean(dat_quad_MM$Burden_sibling_EG_inhRD)
cohensD(dat_quad_MM$Burden_proband_EG_inhRD,dat_quad_MM$Burden_sibling_EG_inhRD,method = "paired") 
wilcox.test(dat_quad_MM$Burden_proband_EG_inhRD, dat_quad_MM$Burden_sibling_EG_inhRD, paired=TRUE, alternative="greater") 

dim(dat_quad_FF) #122
mean(dat_quad_FF$Burden_proband_EG_inhRD) 
mean(dat_quad_FF$Burden_sibling_EG_inhRD) 
cohensD(dat_quad_FF$Burden_proband_EG_inhRD,dat_quad_FF$Burden_sibling_EG_inhRD,method = "paired") 
wilcox.test(dat_quad_FF$Burden_proband_EG_inhRD, dat_quad_FF$Burden_sibling_EG_inhRD, paired=TRUE, alternative="greater") 

## Burden difference between pairs of siblings: inhRD in NEG
dim(dat_quad) #1781
mean(dat_quad$Burden_proband_NEG_inhRD)
mean(dat_quad$Burden_sibling_NEG_inhRD) 
cohensD(dat_quad$Burden_proband_NEG_inhRD,dat_quad$Burden_sibling_NEG_inhRD,method = "paired") 
wilcox.test(dat_quad$Burden_proband_NEG_inhRD, dat_quad$Burden_sibling_NEG_inhRD, paired=TRUE, alternative="greater") 

dim(dat_quad_FM) #101
mean(dat_quad_FM$Burden_proband_NEG_inhRD) 
mean(dat_quad_FM$Burden_sibling_NEG_inhRD) 
cohensD(dat_quad_FM$Burden_proband_NEG_inhRD,dat_quad_FM$Burden_sibling_NEG_inhRD,method = "paired") 
wilcox.test(dat_quad_FM$Burden_proband_NEG_inhRD, dat_quad_FM$Burden_sibling_NEG_inhRD, paired=TRUE, alternative="greater") 

dim(dat_quad_MF) #826
mean(dat_quad_MF$Burden_proband_NEG_inhRD)
mean(dat_quad_MF$Burden_sibling_NEG_inhRD) 
cohensD(dat_quad_MF$Burden_proband_NEG_inhRD,dat_quad_MF$Burden_sibling_NEG_inhRD,method = "paired") 
wilcox.test(dat_quad_MF$Burden_proband_NEG_inhRD, dat_quad_MF$Burden_sibling_NEG_inhRD, paired=TRUE, alternative="greater") 

dim(dat_quad_MM) #732
mean(dat_quad_MM$Burden_proband_NEG_inhRD) 
mean(dat_quad_MM$Burden_sibling_NEG_inhRD) 
cohensD(dat_quad_MM$Burden_proband_NEG_inhRD,dat_quad_MM$Burden_sibling_NEG_inhRD,method = "paired")
wilcox.test(dat_quad_MM$Burden_proband_NEG_inhRD, dat_quad_MM$Burden_sibling_NEG_inhRD, paired=TRUE, alternative="greater") 

dim(dat_quad_FF) #122
mean(dat_quad_FF$Burden_proband_NEG_inhRD) 
mean(dat_quad_FF$Burden_sibling_NEG_inhRD) 
cohensD(dat_quad_FF$Burden_proband_NEG_inhRD,dat_quad_FF$Burden_sibling_NEG_inhRD,method = "paired")
wilcox.test(dat_quad_FF$Burden_proband_NEG_inhRD, dat_quad_FF$Burden_sibling_NEG_inhRD, paired=TRUE, alternative="greater") 

###Table 1.
dat$Burden_proband_EG_all=dat$Burden_proband_EG_dnLoF+dat$Burden_proband_EG_dnNSD+dat$Burden_proband_EG_inhRD
dat$Burden_proband_NEG_all=dat$Burden_proband_NEG_dnLoF+dat$Burden_proband_NEG_dnNSD+dat$Burden_proband_NEG_inhRD
dat_male=dat[which(dat$probandGender=='M'),] #2043 male probands
dat_female=dat[which(dat$probandGender=='F'),] #325 female probands

summary(m1 <- glm(Proband.SRS.Total ~ Burden_proband_EG_all, family="poisson", data=dat_male)) #0.0018599  0.0003814   4.876 1.08e-06 ***
summary(m1 <- glm(Proband.SRS.Total ~ Burden_proband_NEG_all, family="poisson", data=dat_male)) #0.0004072  0.0003244   1.255    0.209
summary(m1 <- glm(Proband.SRS.Total ~ Burden_proband_EG_all, family="poisson", data=dat_female)) #-0.0015107  0.0008772  -1.722    0.085 .
summary(m1 <- glm(Proband.SRS.Total ~ Burden_proband_NEG_all, family="poisson", data=dat_female)) # -0.0030838  0.0006815  -4.525 6.04e-06 ***

png(/home/jixiao/pgfi_data/ifs/jixiao/programs/autismAnalysis/autism_Sanders2015/masterPhenoF/)
ggplot2()


###Table S4.
summary(m1 <- glm(Proband.Verbal.IQ ~ Burden_proband_EG_all, family="poisson", data=dat)) #  -0.007279   0.000400   -18.2   <2e-16 ***
summary(m1 <- glm(Proband.Nonverbal.IQ ~ Burden_proband_EG_all, family="poisson", data=dat)) # -0.0053074  0.0003827  -13.87   <2e-16 ***
summary(m1 <- glm(Proband.Verbal.IQ ~ Burden_proband_NEG_all, family="poisson", data=dat)) # -0.0071722  0.0003361  -21.34   <2e-16 ***
summary(m1 <- glm(Proband.Nonverbal.IQ ~ Burden_proband_NEG_all, family="poisson", data=dat)) # -0.0049064  0.0003202  -15.32   <2e-16 ***
