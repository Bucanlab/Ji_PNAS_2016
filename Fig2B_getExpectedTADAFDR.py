import numpy as np
import random

def getExpectedTADAFDR(TADAf,nPerm,glF): #TADAf should be sorted by TADA FDR value!
	glSet={}
	f0=open(glF)
	f0.readline()
	for l in f0:
		ll=l.rstrip().split('\t')
		if len(ll)>1:
			glSet[ll[0]]=ll[1]
		else:
			glSet[ll[0]]=None
	f0.close()
	glName=glF.split('/')[-1].split('.')[0]
	
	
	outGeneL=[]
	outFDRl=[]
	glNameL=[]
	outF=TADAf+'.'+str(glName)+'_'+str(nPerm)
	f1=open(TADAf)
	header=f1.readline()
	fdrL=[]
	for l in f1:
		ll=l.rstrip().split('\t')
		gene=ll[0]
		value=float(ll[15])
		fdrL.append(value)
		if gene in glSet:
			outGeneL.append(gene)
			outFDRl.append(ll[15])
			if glSet[gene]!=None:
				glNameL.append(glSet[gene])
			else:
				glNameL.append(glName)
	f1.close()
	n=len(fdrL)
	
	nGeneSet=len(outGeneL)
	permMat=[]
	for i in range(nPerm):
		print "Permutation:",i
		permL=sample_wr(fdrL, nGeneSet)
		permL.sort()
		permMat.append(permL)
	
	permMat=np.array(permMat)
	expectedL=np.median(permMat,axis=0)
	expectedL_low=np.percentile(permMat,2.5,axis=0)
	expectedL_high=np.percentile(permMat,97.5,axis=0)
	
	f2=open(outF,'w')
	f2.write('Gene\tTADA_FDR\tTADA_FDR_expected\tTADA_FDR_expected_2.5\tTADA_FDR_expected_97.5\tglName\n')
	for i in range(nGeneSet):
		
		f2.write(outGeneL[i]+'\t'+outFDRl[i]+'\t'+str(expectedL[i])+'\t'+str(expectedL_low[i])+'\t'+str(expectedL_high[i])+'\t'+glNameL[i]+'\n')
	f2.close()
	
def sample_wr(population, k):
    "Chooses k random elements (with replacement) from a population"
    n = len(population)
    _random, _int = random.random, int  # speed hack 
    result = [None] * k
    for i in xrange(k):
        j = _int(_random() * n)
        result[i] = population[j]
    return result	

	
nPerm=100000
TADAf=".../Sanders_2015_S6_TADA.txt"
glF=".../NEG_all_v2.txt"
getExpectedTADAFDR(TADAf,nPerm,glF)
glF=".../EG_all_cell_v2.txt"
getExpectedTADAFDR(TADAf,nPerm,glF)
