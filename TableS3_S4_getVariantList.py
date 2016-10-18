from string import join

def getVariantList_denovo(denovoF1,denovoF2,EGf,NEGf,outF):
	
	EGset=getGeneSet(EGf)
	NEGset=getGeneSet(NEGf)
	f0=open(denovoF1)
	header=f0.readline()
	headerL=header.rstrip().split('\t')
	
	varTypeDic1={"LOF":set(["frame-shift","noEnd","noStart","nonsense","splice-site","no-frame-shift-newStop"]),"NS":set(["frame-shift","noEnd","noStart","nonsense","splice-site","no-frame-shift-newStop","missense"]),"S":set(["frame-shift","noEnd","noStart","nonsense","splice-site","no-frame-shift-newStop","missense","no-frame-shift","synonymous"]),"A":set(["frame-shift","noEnd","noStart","nonsense","splice-site","no-frame-shift-newStop","missense","no-frame-shift","synonymous","intron","3'UTR","non-coding","5'UTR","5'UTR-intron","non-coding-intron","3'UTR-intron"]),"NSonly":set(["missense"])}
	f2=open(outF,'w')
	f2.write("fam\tinChild\tvarType\teffectType\tgeneSet\tsource\t"+'\t'.join(headerL[16:])+'\n') #write header
	for l in f0:
		ll=l.rstrip().split('\t')
		fam=ll[0]
		inChild=ll[4]
		effectGene=ll[6]
		effectType=ll[7]
		
		source="Iossifov et al."
		#familyId,location,variant,vcfVariant,inChild,fromParent,effectGene,effectType,familyDescription,CSHL,YALE,UW,IossifovWE2012,EichlerWE2012,StateWE2012,EichlerTG2012=ll[0:]
		if effectGene in EGset:
			geneSet="EG"
		elif effectGene in NEGset:
			geneSet="NEG"
		else:
			geneSet="Unknown"
		if effectType in varTypeDic1["LOF"]:
			varType="dnLoF"
		elif effectType in varTypeDic1["NSonly"]:
			varType="dnNSD"
		else:
			continue
			
		if varType=="dnNSD":
			cadd=ll[66]
			if cadd=='.':
				cadd_indel=ll[68]
				if cadd_indel=='.':
					continue
				elif float(cadd_indel)<10:
					continue
				else:
					pass		
			elif float(cadd)<10:
				continue
			else:
				pass		
		outL=[fam,inChild,varType,effectType,geneSet,source]+ll[16:]
		f2.write('\t'.join(outL)+'\n')
		
	f0.close()
	
	f1=open(denovoF2)
	header=f1.readline()
	headerL=header.rstrip().split('\t')
	varTypeDic2={"LOF":set(["STOP_GAINED","FRAME_SHIFT","SPLICE_SITE_ACCEPTOR","START_LOST","STOP_LOST","SPLICE_SITE_DONOR"]),"NS":set(["STOP_GAINED","FRAME_SHIFT","SPLICE_SITE_ACCEPTOR","START_LOST","STOP_LOST","SPLICE_SITE_DONOR","NON_SYNONYMOUS_CODING","NON_SYNONYMOUS_START","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION"]),"S":set(["STOP_GAINED","FRAME_SHIFT","SPLICE_SITE_ACCEPTOR","START_LOST","STOP_LOST","SPLICE_SITE_DONOR","NON_SYNONYMOUS_CODING","NON_SYNONYMOUS_START","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","SYNONYMOUS_CODING","SYNONYMOUS_STOP","EXON"]),"A":set(["STOP_GAINED","FRAME_SHIFT","SPLICE_SITE_ACCEPTOR","START_LOST","STOP_LOST","SPLICE_SITE_DONOR","NON_SYNONYMOUS_CODING","NON_SYNONYMOUS_START","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","SYNONYMOUS_CODING","SYNONYMOUS_STOP","EXON","INTRON","UTR_3_PRIME","UTR_5_PRIME"]),"NSonly":set(["NON_SYNONYMOUS_CODING","NON_SYNONYMOUS_START","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION"])}
	
	for l in f1:
		ll=l.rstrip().split('\t')
		inChild=ll[4].split('.')[-1]
		effectGene=ll[5]
		effectType=ll[6]
		fam=ll[4].split('.')[0]
		source="Krumm et al."
		if effectGene in EGset:
			geneSet="EG"
		elif effectGene in NEGset:
			geneSet="NEG"
		else:
			geneSet="Unknown"
		if effectType in varTypeDic2["LOF"]:
			varType="dnLoF"
		elif effectType in varTypeDic2["NSonly"]:
			varType="dnNSD"
		else:
			continue
		if varType=="dnNSD":
			cadd=ll[78]
			if cadd=='.':
				cadd_indel=ll[80]
				if cadd_indel=='.':
					continue
				elif float(cadd_indel)<10:
					continue
				else:
					pass
			elif float(cadd)<10:
				continue
			else:
				pass
		outL=[fam,inChild,varType,effectType,geneSet,source]+ll[28:]
		f2.write('\t'.join(outL)+'\n')
		
	f1.close()
	f2.close()	

def getVariantList_inherited(inheritedF,EGf,NEGf,outF):
	EGset=getGeneSet(EGf)
	NEGset=getGeneSet(NEGf)
	f1=open(inheritedF)
	header=f1.readline()
	headerL=header.rstrip().split('\t')
	f2=open(outF,'w')
	f2.write("fam\tinChild\tvarType\teffectType\tgeneSet\tsource\t"+'\t'.join(headerL[0:18])+'\t'+headerL[70]+'\t'+headerL[88]+'\n')
	
	varTypeDic2={"LOF":set(["STOP_GAINED","FRAME_SHIFT","SPLICE_SITE_ACCEPTOR","START_LOST","STOP_LOST","SPLICE_SITE_DONOR"]),"NS":set(["STOP_GAINED","FRAME_SHIFT","SPLICE_SITE_ACCEPTOR","START_LOST","STOP_LOST","SPLICE_SITE_DONOR","NON_SYNONYMOUS_CODING","NON_SYNONYMOUS_START","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION"]),"S":set(["STOP_GAINED","FRAME_SHIFT","SPLICE_SITE_ACCEPTOR","START_LOST","STOP_LOST","SPLICE_SITE_DONOR","NON_SYNONYMOUS_CODING","NON_SYNONYMOUS_START","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","SYNONYMOUS_CODING","SYNONYMOUS_STOP","EXON"]),"A":set(["STOP_GAINED","FRAME_SHIFT","SPLICE_SITE_ACCEPTOR","START_LOST","STOP_LOST","SPLICE_SITE_DONOR","NON_SYNONYMOUS_CODING","NON_SYNONYMOUS_START","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","SYNONYMOUS_CODING","SYNONYMOUS_STOP","EXON","INTRON","UTR_3_PRIME","UTR_5_PRIME"]),"NSonly":set(["NON_SYNONYMOUS_CODING","NON_SYNONYMOUS_START","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION"])}
	
	for l in f1:
		ll=l.rstrip().split('\t')
		inChild=ll[9]
		effectGene=ll[11]
		effectType=ll[12]
		fam=ll[6]
		Set=ll[10]
		source="Krumm et al."
		if Set!="Intersection": # Use variants called by both FreeBayes and GATK
			continue
		
		if effectGene in EGset:
			geneSet="EG"
		elif effectGene in NEGset:
			geneSet="NEG"
		elif effectGene !="":
			geneSet="Unknown"
		else:
			continue
		if effectType in varTypeDic2["LOF"]:
			varType="inhRD"
		elif effectType in varTypeDic2["NSonly"]:
			varType="inhRD"
		else:
			continue
		if varType=="inhRD":
			cadd=ll[88]
			if cadd=='.':
				continue
			elif float(cadd)<10:
				continue
			else:
				pass
			MAF=ll[70]
			if MAF=='':
				pass
			elif float(MAF)<=0.01:
				pass
			else:
				continue
		else:
			continue	
		outL=[fam,inChild,varType,effectType,geneSet,source]+ll[0:18]+[ll[70]]+[ll[88]]
		f2.write('\t'.join(outL)+'\n')
		
	f1.close()
	f2.close()	

def addMutationBurden(masterPhenoF,inheritF,denovoF,EGf,NEGf,outF):
	
	f0=open(masterPhenoF)
	f1=open(outF,'w')
	header=f0.readline()
	f1.write(header.rstrip()+'\tBurden_proband_EG_dnLoF\tBurden_sibling_EG_dnLoF\tBurden_proband_EG_dnNSD\tBurden_sibling_EG_dnNSD\tBurden_proband_EG_inhRD\tBurden_sibling_EG_inhRD\tBurden_proband_NEG_dnLoF\tBurden_sibling_NEG_dnLoF\tBurden_proband_NEG_dnNSD\tBurden_sibling_NEG_dnNSD\tBurden_proband_NEG_inhRD\tBurden_sibling_NEG_inhRD\n')
	
	EGset=getGeneSet(EGf)
	NEGset=getGeneSet(NEGf)
	
	famDic={}
	f2=open(inheritF)
	f2.readline()
	for l in f2:
		ll=l.rstrip().split('\t')
		famDic.setdefault(ll[0],[0]*12)
		gene=ll[17]
		inChild=ll[1]
		if gene in EGset:
			if inChild.endswith("both"):
				famDic[ll[0]][4]+=1
				famDic[ll[0]][5]+=1
			elif inChild.endswith("pro"):
				famDic[ll[0]][4]+=1
			elif inChild.endswith("sib"):
				famDic[ll[0]][5]+=1
		elif gene in NEGset:
			if inChild.endswith("both"):
				famDic[ll[0]][10]+=1
				famDic[ll[0]][11]+=1
			elif inChild.endswith("pro"):
				famDic[ll[0]][10]+=1
			elif inChild.endswith("sib"):
				famDic[ll[0]][11]+=1
	f2.close()		
	
	f3=open(denovoF)
	f3.readline()
	for l in f3:
		ll=l.rstrip().split('\t')
		famDic.setdefault(ll[0],[0]*12)
		gene=ll[12]
		inChild=ll[1]
		varType=ll[2]
		if varType=="dnLoF":
			if gene in EGset:
				if len(inChild)==4:
					famDic[ll[0]][0]+=1
					famDic[ll[0]][1]+=1
				elif inChild.startswith("p"):
					famDic[ll[0]][0]+=1
				elif inChild.startswith("s"):
					famDic[ll[0]][1]+=1
			elif gene in NEGset:
				if len(inChild)==4:
					famDic[ll[0]][6]+=1
					famDic[ll[0]][7]+=1
				elif inChild.startswith("p"):
					famDic[ll[0]][6]+=1
				elif inChild.startswith("s"):
					famDic[ll[0]][7]+=1
		elif varType=="dnNSD":
			if gene in EGset:
				if len(inChild)==4:
					famDic[ll[0]][2]+=1
					famDic[ll[0]][3]+=1
				elif inChild.startswith("p"):
					famDic[ll[0]][2]+=1
				elif inChild.startswith("s"):
					famDic[ll[0]][3]+=1
			elif gene in NEGset:
				if len(inChild)==4:
					famDic[ll[0]][8]+=1
					famDic[ll[0]][9]+=1
				elif inChild.startswith("p"):
					famDic[ll[0]][8]+=1
				elif inChild.startswith("s"):
					famDic[ll[0]][9]+=1
	f3.close()
	
	for l in f0:
		ll=l.rstrip().split('\t')
		if ll[5]=="trio":
			for i in [1,3,5,7,9,11]:
				famDic[ll[0]][i]=""
		outL=ll+famDic[ll[0]]
		outL=[str(x) for x in outL]
		f1.write('\t'.join(outL)+'\n')
	f0.close()
	f1.close()	
	

#Recent protocols
denovoF1="/home/jixiao/pgfi_data/ifs/jixiao/programs/autismAnalysis/autism_Sanders2015/annovar/Iossifov_2014_denovo.hg19_multianno.txt2"
denovoF2="/home/jixiao/pgfi_data/ifs/jixiao/programs/autismAnalysis/autism_Sanders2015/annovar/Krumm_2015_denovo.hg19_multianno.txt2"
inheritedF="/home/jixiao/pgfi_data/ifs/jixiao/programs/autismAnalysis/autism_Sanders2015/Krumm_2015_All_Rare_Inherited_Final.txt"
EGf="/home/jixiao/pgfi_data/ifs/jixiao/programs/MONSTER/geneLists/EG_all_cell_v2.txt"
NEGf="/home/jixiao/pgfi_data/ifs/jixiao/programs/MONSTER/geneLists/NEG_all_v2.txt"
outF="/home/jixiao/pgfi_data/ifs/jixiao/programs/autismAnalysis/autism_Sanders2015/varList_dnLoF_dnNSD.txt"
#getVariantList_denovo(denovoF1,denovoF2,EGf,NEGf,outF)
outF="/home/jixiao/pgfi_data/ifs/jixiao/programs/autismAnalysis/autism_Sanders2015/varList_inhRD.txt"
#getVariantList_inherited(inheritedF,EGf,NEGf,outF)

masterPhenoF="/home/jixiao/pgfi_data/ifs/jixiao/programs/autismAnalysis/autism_Sanders2015/masterPhenoF/famPhenoMasterTable_shared_curated.txt"
outF="/home/jixiao/pgfi_data/ifs/jixiao/programs/autismAnalysis/autism_Sanders2015/masterPhenoF/masterPheno_indivBurden_final.txt"
inheritF="/home/jixiao/pgfi_data/ifs/jixiao/programs/autismAnalysis/autism_Sanders2015/varList_inhRD.txt"
denovoF="/home/jixiao/pgfi_data/ifs/jixiao/programs/autismAnalysis/autism_Sanders2015/varList_dnLoF_dnNSD.txt"
EGf="/home/jixiao/pgfi_data/ifs/jixiao/programs/MONSTER/geneLists/EG_all_cell_v2.txt"
NEGf="/home/jixiao/pgfi_data/ifs/jixiao/programs/MONSTER/geneLists/NEG_all_v2.txt"
#addMutationBurden(masterPhenoF,inheritF,denovoF,EGf,NEGf,outF)

masterPhenoF="/home/jixiao/pgfi_data/ifs/jixiao/programs/autismAnalysis/autism_Sanders2015/masterPhenoF/famPhenoMasterTable_shared_curated.txt"
outF="/home/jixiao/pgfi_data/ifs/jixiao/programs/autismAnalysis/autism_Sanders2015/masterPhenoF/masterPheno_indivBurden_noSFARI.txt"
inheritF="/home/jixiao/pgfi_data/ifs/jixiao/programs/autismAnalysis/autism_Sanders2015/varList_inhRD.txt"
denovoF="/home/jixiao/pgfi_data/ifs/jixiao/programs/autismAnalysis/autism_Sanders2015/varList_dnLoF_dnNSD.txt"
EGf="/home/jixiao/pgfi_data/ifs/jixiao/programs/MONSTER/geneLists/EG_all_cell_v2_noSFARI.txt"
NEGf="/home/jixiao/pgfi_data/ifs/jixiao/programs/MONSTER/geneLists/NEG_all_v2_noSFARI.txt"
addMutationBurden(masterPhenoF,inheritF,denovoF,EGf,NEGf,outF)
