#!/usr/bin/env python
# coding: latin1

import random
import math
import os
import sys
import time
import numpy as np

# Simulated data assessment: 

# Assumes that the working directory contains the files of one data set.
# And the calculated clusterings are in another file passed as parameter.


# First we load all data:
cwd = os.getcwd()

########## The copy numbers of the reads:
readcopynumberpath=""
for path in os.listdir(cwd):
	print path
	if path[-len("_ReadCopynumbers"):]=="_ReadCopynumbers":
		readcopynumberpath=path
if readcopynumberpath=="":
	print "readcopynumberpath not detected."
	sys.exit()
if os.path.exists(readcopynumberpath):
	Read2Copy=[]
	f=open(readcopynumberpath,'r')
	for line in f.readlines():
		Read2Copy.append(int(line))
	f.close()
else:
	print "readcopynumberpath missing."
	sys.exit()
print "len(Read2Copy): {}".format(len(Read2Copy))
#######################################


################# Ground truth placement of reads
readplacementpath=""
for path in os.listdir(cwd):
	if path[-len("_ReadPlacements"):]=="_ReadPlacements":
		readplacement=path
if readplacement=="":
	print "readplacement not detected."
	sys.exit()
if os.path.exists(readplacement):
	Read2Place=[]
	f=open(readplacement,'r')
	for line in f.readlines():
		Read2Place.append(int(line))
print "len(Read2Place): {}".format(len(Read2Place))
################################################


###################### ReadSeqInfo
readseqinfopath=""
for path in os.listdir(cwd):
	if path[-len("_ReadSeqInfo"):]=="_ReadSeqInfo":
		readseqinfopath=path
if readseqinfopath=="":
	print "readseqinfopath not detected."
	sys.exit()

if os.path.exists(readseqinfopath):
	Seq2Reads=[]
	f=open(readseqinfopath,'r')
	read=0
	for line in f.readlines():
		liste=line.split()
		for x in liste:
			Seq2Reads.append(read)
		read+=1
	f.close()
else:
	print "readseqinfopath error."
	sys.exit()

print "len(Seq2Reads): {}".format(len(Seq2Reads))

######################################


################### The sequences ######
#"_Seq.fasta"  # Sequences
seqpath=""
for path in os.listdir(cwd):
	if path[-len("_Seq.fasta"):]=="_Seq.fasta":
		seqpath=path
if seqpath=="":
	print "seqpath not detected."
	sys.exit()

if os.path.exists(seqpath):	
	f=open(seqpath,'r')
	Seqs=[]
	read=''
	for line in f.readlines():
		if line[0]=='>' and read!='':
			Seqs.append(read)
			read=''
		else:
			read+=line[:-1]
	Seqs.append(read)
	f.close()

else:
	print "seqpath error."
	sys.exit()
########################################


#################### The Reads #########
".fasta"  # Reads
readpath=""
for path in os.listdir(cwd):
	if path[-len(".fasta"):]==".fasta" and path[-len("_Seq.fasta"):]!="_Seq.fasta" and path[-len("_Template.fasta"):]!="_Template.fasta":
		readpath=path
if readpath=="":
	print "readpath not detected."
	sys.exit()

if os.path.exists(readpath):
	f=open(readpath,'r')
	Reads=[]
	read=''
	for line in f.readlines():
		if line[0]=='>' and read!='':
			Reads.append(read)
			read=''
		else:
			read+=line[:-1]
	Reads.append(read)
	f.close()
else:
	print "readpath error."
	sys.exit()
#####################################################


##################### The seq classes, MSA and unique sequence.
seqclasspath=""
datasetname=""
for path in os.listdir(cwd):
	if path[-len("_SeqClass"):]=="_SeqClass":
		seqclasspath=path
		datasetname=path[:-len("_SeqClass")]
if seqclasspath=="":
	print "seqclasspath not detected."
	sys.exit()

if os.path.exists(seqclasspath):
	MSA2Seq=[]
	Unique2Seq=[]
	f=open(seqclasspath,'r')
	seq=0
	for line in f.readlines():
		if line[:1]=='r':
			MSA2Seq.append(seq)
		else:
			Unique2Seq.append(seq)
		seq+=1
	f.close()
else:
	print "seqclasspath error."
	sys.exit()
if datasetname=="":
	print "datasetname could not be extracted from seqclasspath."
	sys.exit()
print "len(MSA2Seq): {}".format(len(MSA2Seq))
print "len(Unique2Seq): {}".format(len(Unique2Seq))
#############################################



#################### The resolution path
if len(sys.argv)<2:
	print "Missing resolution path."
	sys.exit()

resolutionpath=sys.argv[1]

if os.path.exists(resolutionpath):
	Resolutions=[]
	for path in os.listdir(resolutionpath):
		if path[:len("KmeansSubdivisionOf_")]=="KmeansSubdivisionOf_":
			liste=path.split('_')
			start=int(liste[1])
			Resolution=[]
			print resolutionpath+path
			f=open(resolutionpath+path,'r')
			for line in f.readlines():
				Resolution.append(int(line))
			f.close()
			Resolutions.append((start,Resolution))
	print "Resolution starts:",
	print [no for (no,Res) in sorted(Resolutions)]
	Resolutions=[Res for (no,Res) in sorted(Resolutions)]
else:
	print "resolutionpath error"
	sys.exit()
print "Lengths of Resolutions:"
for res in Resolutions:
	print len(res),
print
#########################################

##### Here we calculate the ReadSeqInfo again:
Seq2Reads=[]

r=0
s=0
def SeqInRead(seq,read):
	return seq[:100] in read

while s<len(Seqs) and r<len(Reads):
	if SeqInRead(Seqs[s],Reads[r]):
		s+=1
		Seq2Reads.append(r)
	else:
		r+=1

############# Calculation of flanking clusters:
FlankingRight=[]
for t in range(len(MSA2Seq)): # Durch die signaturen, rows, ...
	if MSA2Seq[t]+1 in Unique2Seq and Seq2Reads[MSA2Seq[t]+1]==Seq2Reads[MSA2Seq[t]]:
		FlankingRight.append(Read2Copy[Seq2Reads[MSA2Seq[t]]])
	else:
		FlankingRight.append(-1)
print "MSA2Seq[-1]",
print MSA2Seq[-1]
FlankingLeft=[]
for t in range(len(MSA2Seq)): # Durch die signaturen, rows, ...
	if MSA2Seq[t]-1 in Unique2Seq and Seq2Reads[MSA2Seq[t]-1]==Seq2Reads[MSA2Seq[t]]:
		FlankingLeft.append(Read2Copy[Seq2Reads[MSA2Seq[t]]])
	else:
		FlankingLeft.append(-1)

print "{}/{} FlankingLeft".format(sum(FlankingLeft),len(FlankingLeft))
print "{}/{} FlankingRight".format(sum(FlankingRight),len(FlankingRight))

print "Lengths of Flankings: {} {}".format(len(FlankingLeft),len(FlankingRight))
print "datasetname: {}".format(datasetname)

# f=open(datasetname+'_FlankingLeft','w')
# for x in FlankingLeft:
# 	f.write(str(x)+'\n')
# f.close()

# f=open(datasetname+'_FlankingRight','w')
# for x in FlankingRight:
# 	f.write(str(x)+'\n')
# f.close()

############################################################


# Then we assess single resolutions:

# First a ground truth for alle MSASeq: In three steps: MSA2Seq+Seq2Reads+Read2Copy
GroundTruthResolution=[Read2Copy[Seq2Reads[MSA2Seq[z]]] for z in range(len(MSA2Seq))]

def GroupMaker(Resolution):
	return [[z for z in range(len(Resolution)) if Resolution[z]==x] for x in range(max(Resolution)+1) if Resolution.count(x)>0]

def ResolutionQuality(GroundTruthResolution, Resolution):
	FlankingLeft=[-1 for z in range(len(Resolution))]
	FlankingRight=[-1 for z in range(len(Resolution))]
	GroundTruthResolution2=[]  # Only those that occur in the resolution under assessment
	for z in range(len(GroundTruthResolution)):
		if Resolution[z]>-1:
			GroundTruthResolution2.append(GroundTruthResolution[z])
		else:
			GroundTruthResolution2.append(-1) 

	GroundTruthGroups=GroupMaker(GroundTruthResolution2)
	for g in range(len(GroundTruthGroups)):
		for z in GroundTruthGroups[g]:
			FlankingLeft[z]=g
			FlankingRight[z]=g

	Matrix1=[[0.0 for ttt in range(max(Resolution)+1)] for tt in range(len(GroundTruthGroups))]
	Matrix2=[[0.0 for ttt in range(len(GroundTruthGroups))] for tt in range(max(Resolution)+1)]

	Matrix1=np.array(Matrix1)
	Matrix2=np.array(Matrix2)

	for tt in range(len(GroundTruthGroups)):
		size=float(len(GroundTruthGroups[tt]))
		for ttt in range(max(Resolution)+1): 
			Matrix1[tt][ttt]=float(len([z for z in GroundTruthGroups[tt] if Resolution[z]==ttt]))/size

	for ttt in range(max(Resolution)+1):
		size=float(Resolution.count(ttt))
		for tt in range(len(GroundTruthGroups)):
			if size>0:
				Matrix2[ttt][tt]=float(len([z for z in GroundTruthGroups[tt] if Resolution[z]==ttt]))/size

	Matrix3=np.dot(Matrix1,Matrix2)

	# Normalisation: 
	for tt in range(len(Matrix3)):
		summe=sum([Matrix3[tt][ttt] for ttt in range(len(Matrix3[tt]))])
		for ttt in range(len(Matrix3[tt])):
			if summe>0.0:
				Matrix3[tt][ttt]/=summe

	# for conconf in [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
	# 	count=0
	# 	for tt in range(len(Matrix3)):
	# 		if Matrix3[tt][tt]==max(Matrix3[tt]) and max(Matrix3[tt])>conconf:
	# 			count+=1

	#  	print "Bei conconf {}: {} falsely resolved.".format(float(conconf)/10,sum([1 for tt in range(len(Matrix3)) if (Matrix3[tt][tt]!=max(Matrix3[tt]) and max(Matrix3[tt])>float(conconf)/10) and max(Matrix3[ np.nonzero(Matrix3[tt] == max(Matrix3[tt]))[0][0] ])==max(Matrix3[tt])]))

	# #print "Alternative Berechnung:"
	conconfpositives=[0 for c in range(10)]
	truepositives=0
	falsepositives=0
	theresolved=[]
	maxis=[]
	for t in range(len(Matrix3)):
		maxi=0.0
		maxtt=0
		for tt in range(len(Matrix3)):
			if Matrix3[t][tt]>maxi:
				maxi=Matrix3[t][tt]
				maxtt=tt
		if maxi==max(Matrix3[maxtt]):
			if maxtt!=t:
				falsepositives+=1
				theresolved.append(-1)
			if maxtt==t:
				theresolved.append(1)
				truepositives+=1
				for c in range(10):  # true positives + >conconfcutoff
					if maxi>float(c)/10.0:
						conconfpositives[c]+=1
		else:
			theresolved.append(0)
		maxis.append(maxi)

	print "truepositives {}, falsepositives {} bei cutoff 0.0.".format(truepositives,falsepositives)
	print "Number of resolved copies by cutoff > 0.0, 0.1, 0.2 ... 0.9:"
	print conconfpositives
	print 

	return theresolved,maxis


for Resolution in Resolutions:
	ResolutionQuality(GroundTruthResolution, Resolution)

# Then we assess multi-step resolutions:

def ProbabilityMatrix(Resolution1,Resolution2):
	Matrix=[[0.0 for t in range(max(Resolution2)+1)] for tt in range(max(Resolution1)+1)]
	Sums=[len([tt for tt in range(len(Resolution1)) if Resolution1[tt]==t and Resolution2[tt]>-1]) for t in range(max(Resolution1)+1)]
	for t in range(len(Resolution1)):
		if Resolution1[t]>-1 and Resolution2[t]>-1:
			Matrix[Resolution1[t]][Resolution2[t]]+=1.0
	for t in range(max(Resolution1)+1):
		for tt in range(max(Resolution2)+1):
			if Sums[t]>0:
				Matrix[t][tt]/=float(Sums[t])
	Matrix=np.array(Matrix)
	return Matrix

def MultiStepResolution(FlankingLeft,Resolutions,FlankingRight):
	ForwardMatrices=[]
	BackwardMatrices=[]
	AllResolutions=[FlankingLeft]+Resolutions+[FlankingRight]
	for r in range(len(AllResolutions)-1):
		ForwardMatrices.append(ProbabilityMatrix(AllResolutions[r],AllResolutions[r+1]))
		BackwardMatrices.append(ProbabilityMatrix(AllResolutions[len(AllResolutions)-1-r],AllResolutions[len(AllResolutions)-2-r]))
	ForwardConCon=np.dot(ForwardMatrices[0],ForwardMatrices[1])
	BackwardConCon=np.dot(BackwardMatrices[0],BackwardMatrices[1])
	for t in range(2,len(ForwardMatrices)):
		ForwardConCon=np.dot(ForwardConCon,ForwardMatrices[t])
		BackwardConCon=np.dot(BackwardConCon,BackwardMatrices[t])
	AllConCon=np.multiply(ForwardConCon,np.transpose(BackwardConCon))

	# Normalisation: 
	for tt in range(len(AllConCon)):
		summe=sum([AllConCon[tt][ttt] for ttt in range(len(AllConCon[tt]))])
		for ttt in range(len(AllConCon[tt])):
			if summe>0.0:
				AllConCon[tt][ttt]/=summe

	# #print "Alternative Berechnung:"
	conconfpositives=[0 for c in range(10)]
	truepositives=0
	falsepositives=0
	theresolved=[]
	maxis=[]
	for t in range(len(AllConCon)):
		maxi=0.0
		maxtt=0
		for tt in range(len(AllConCon)):
			if AllConCon[t][tt]>maxi:
				maxi=AllConCon[t][tt]
				maxtt=tt
		if maxi==max(AllConCon[maxtt]):
			if maxtt!=t:
				falsepositives+=1
				theresolved.append(-1)
			if maxtt==t:
				theresolved.append(1)
				truepositives+=1
				for c in range(10):  # true positives + >conconfcutoff
					if maxi>float(c)/10.0:
						conconfpositives[c]+=1
		else:
			theresolved.append(0)
		maxis.append(maxi)

	print "MultiStepResolution:"
	print "truepositives {}, falsepositives {} bei cutoff 0.0.".format(truepositives,falsepositives)
	print "Number of resolved copies by cutoff > 0.0, 0.1, 0.2 ... 0.9:"
	print conconfpositives
	print 

	return AllConCon


MultiStepResolution(FlankingLeft,Resolutions,FlankingRight)




sys.exit()









