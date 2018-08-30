#!/usr/bin/env python
# coding: latin1

import random
import math
import os
import sys
import time
import numpy as np

# This is a tool that takes in a TransposonClustering 
# uses the pathname to load other relevant data
# and calculates interesting stats on the quality of the resolution.



if len(sys.argv)<2:
	print "Missing KmeansSubdivision path."
	sys.exit()

kmeanspath=sys.argv[1]
if not os.path.exists(kmeanspath):
	print "Missing KmeansSubdivision file."
	sys.exit()

liste=kmeanspath.split('_')
start=int(liste[1])
ende=int(liste[2])
number=int(liste[4])


pathsuffix=kmeanspath[19:]

print "Path suffix: {} -> start {}, ende {}, data set number {}".format(pathsuffix,start,ende,number)
print 

maxcorrpath="MaxCorrsOf_MidTransposonMMA_"+str(number)+"_real"
dropoffpath="DropoffSubdivisionOf"+pathsuffix
reldroppath="RelDropSubdivisionOf"+pathsuffix
msapath="MidTransposonMMA_"+str(number)+"_real"
groundtruthpath="TransposonCopies_"+str(number)

Paths=[kmeanspath,maxcorrpath,dropoffpath,reldroppath,msapath,groundtruthpath]

print "Are all files available: "
for path in Paths:
	print "{} : exists {}".format(path,os.path.exists(kmeanspath))
print


def MaxCorrLoader(path):
	Liste=[]
	if os.path.exists(path):
		f=open(path,'r')
		for line in f.readlines():
			if len(line)>1:		
				Liste.append(float(line))
	return Liste

def MSALoader(path):
	MA=[]
	if os.path.exists(path):
		f=open(path,'r')
		for line in f.readlines():
			if len(line)>1:		
				MA.append(line[:-1])
	return MA

def ResolutionLoader(path):
	Resolution=[]
	if os.path.exists(path):
		f=open(path,'r')
		for line in f.readlines():
			if line=='\n':
				Resolution.append(None)
			else:
				Resolution.append(int(line))
		f.close()
	return Resolution

Basen2Index={'a':0,'c':1,'g':2,'t':3,'-':4,' ':5,'A':0,'C':1,'G':2,'T':3,'_':4}
def Konsensus(liste):
	if liste==[]:
		return ''
	kon=''
	for s in range(len(liste[0])):
		counter=np.array([0,0,0,0,0,0])
		for sig in liste:
			counter[Basen2Index[sig[s]]]+=1
		counter[5]=0
		kon+='acgt- '[np.argmax(counter)]
	return kon

def Diff(sig1,sig2):
	return len([1 for t in range(len(sig1)) if sig1[t]!=sig2[t] and ' '!=sig2[t] and sig1[t]!=' '])

def Resolvability(GroundTruthResolution,Signatures):
	GroundTruthGroups=GroupMaker(GroundTruthResolution)
	Kons=[Konsensus([Signatures[z] for z in group]) for group in GroundTruthGroups if len(group)>0]
	summe=[0,0,0,0,0,0,0,0,0,0,0]
	MinDiffs=[]
	for k in range(len(Kons)):
		unique=[1,1,1,1,1,1,1,1,1,1,1]
		mindiff=1000000
		for kk in range(len(Kons)):
			if k!=kk:
				diff=Diff(Kons[k],Kons[kk])
				if diff<mindiff:
					mindiff=diff
				for t in range(diff,len(unique)):
					unique[t]=0
		MinDiffs.append(diff)
		for t in range(len(unique)):
			summe[t]+=unique[t]
	print "According to the statistically significant differences,"
	print "If we demand > 0,1,2,3,4,5,6,7,8,9 differences between copy groups consensuses:"
	print "This data set has {} / {} unique groups.".format(summe,len(Kons))
	print 
	return MinDiffs

def HalfResolvability(GroundTruthResolution,Signatures):
	GroundTruthGroups=GroupMaker(GroundTruthResolution)
	Kons=[Konsensus([Signatures[z] for z in group]) for group in GroundTruthGroups if len(group)>0]

	biggerbigger=0
	biggersmaller=0
	smallerbigger=0
	smallersmaller=0
	MinDiffs1=[]
	MinDiffs2=[]
	for k in range(len(Kons)):
		mindiff1=10000
		mindiff2=10000
		for kk in range(len(Kons)):
			if k!=kk:
				diff1=Diff(Kons[k][:len(Kons[k])/2],Kons[kk][:len(Kons[k])/2])
				diff2=Diff(Kons[k][len(Kons[k])/2:],Kons[kk][len(Kons[k])/2:])
				if diff1<mindiff1:
					mindiff1=diff1
				if diff2<mindiff2:
					mindiff2=diff2

		if mindiff1>5 and mindiff2>5:
			biggerbigger+=1
		if mindiff1>5 and mindiff2<=5:
			biggersmaller+=1
		if mindiff1<=5 and mindiff2>5:
			smallerbigger+=1
		if mindiff1<=5 and mindiff2<=5:
			smallersmaller+=1
		MinDiffs1.append(mindiff1)
		MinDiffs2.append(mindiff2)
	print "{}/{} bigger and {}/{} smaller 5 have >5 in the second half.".format(biggerbigger,biggerbigger+biggersmaller,smallerbigger,smallerbigger+smallersmaller)
	return MinDiffs1,MinDiffs2

def SignaturesMaker(MA,MaxCorrs,cutoff,start,ende):
	return [''.join([MA[z][x] for x in range(start/5,ende/5) if MaxCorrs[x]>cutoff]) for z in range(len(MA))]

def GroupMaker(Resolution):
	return [[z for z in range(len(Resolution)) if Resolution[z]==x] for x in range(max(Resolution)+1) if Resolution.count(x)>0]

def ResolutionQuality(GroundTruthResolution, Resolution):
	FlankingLeft=[-1 for z in range(len(Resolution))]
	FlankingRight=[-1 for z in range(len(Resolution))]
	GroundTruthGroups=GroupMaker(GroundTruthResolution)
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

	count=0
	for tt in range(len(Matrix3)):
		if Matrix3[tt][tt]==max(Matrix3[tt]) and max(Matrix3[tt])>0.1:
			count+=1

	# print "This many resolved:",
	# print count
	# print "Average connection confidence:",
	# print sum([Matrix3[tt][tt] for tt in range(len(Matrix3))])/float(len(Matrix3))

	# for conconf in range(10):
	# 	print "Bei conconf {}: {} resolved.".format(float(conconf)/10,sum([1 for tt in range(len(Matrix3)) if (Matrix3[tt][tt]==max(Matrix3[tt]) and Matrix3[tt][tt]>float(conconf)/10) ]))

	# print "Resolved if conconf>2*second best:",
	# print sum([1 for tt in range(len(Matrix3)) if Matrix3[tt][tt]==max(Matrix3[tt]) and max(Matrix3[tt])>sorted(Matrix3[tt])[-2]*2.0])

	# print "False positives:"
	# for conconf in range(10):
	# 	print "Bei conconf {}: {} falsely resolved.".format(float(conconf)/10,sum([1 for tt in range(len(Matrix3)) if (Matrix3[tt][tt]!=max(Matrix3[tt]) and max(Matrix3[tt])>float(conconf)/10) ]))
	# print "False positives if conconf>2*second best:",
	# print sum([1 for tt in range(len(Matrix3)) if Matrix3[tt][tt]!=max(Matrix3[tt]) and max(Matrix3[tt])>sorted(Matrix3[tt])[-2]*2.0])

	# print "Real false positives:"
	# for conconf in range(10):
	# 	print "Bei conconf {}: {} falsely resolved.".format(float(conconf)/10,sum([1 for tt in range(len(Matrix3)) if (Matrix3[tt][tt]!=max(Matrix3[tt]) and max(Matrix3[tt])>float(conconf)/10) and max(Matrix3[ np.nonzero(Matrix3[tt] == max(Matrix3[tt]))[0][0] ])==max(Matrix3[tt])]))

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


###########
# Loading data
#

MaxCorrs=MaxCorrLoader(maxcorrpath)
MaxCorrs=[max([MaxCorrs[t+i] for i in range(5)]) for t in range(0,len(MaxCorrs),5)] # Reduction onto columns

KmeansResolution=ResolutionLoader(kmeanspath)
DropoffResolution=ResolutionLoader(dropoffpath)
RelDropResolution=ResolutionLoader(reldroppath)
MSA=MSALoader(msapath)
GroundTruthResolution=ResolutionLoader(groundtruthpath)


###########
# Calculation stats
# 
cutoff=1.0  # Eventuell ändern
Signatures=SignaturesMaker(MSA,MaxCorrs,cutoff,start,ende)

# MinDiffs1,MinDiffs2=HalfResolvability(GroundTruthResolution,Signatures)
# print MinDiffs1
# print MinDiffs2
# #sys.exit()

MinDiffs=Resolvability(GroundTruthResolution,Signatures)


#ResolutionQuality(GroundTruthResolution, DropoffResolution)

#print "RelDropResolution:"
#ResolutionQuality(GroundTruthResolution, RelDropResolution)
print "DropoffResolution:"
theresolved,maxis=ResolutionQuality(GroundTruthResolution, DropoffResolution)

print "RelDropResolution:"
theresolved,maxis=ResolutionQuality(GroundTruthResolution, RelDropResolution)

print "KmeansResolution:"
theresolved,maxis=ResolutionQuality(GroundTruthResolution, KmeansResolution)

# Hier vergleiche ich die MinDiffs mit den theresolved. 
# Pro diff wie viele sind resolved, falsepositive or unresolved + gesamtanzahl
# print MinDiffs
# for diff in range(20):
# 	res=0
# 	fal=0
# 	unr=0
# 	cou=0
# 	for t in range(len(MinDiffs)):
# 		if MinDiffs[t]>diff:
# 			if theresolved[t]==1:
# 				res+=1
# 			if theresolved[t]==-1:
# 				fal+=1
# 			if theresolved[t]==0:
# 				unr+=1
# 			cou+=1
# 	print "If mindiff>{}: {} resolved, {} false pos, {} unresolved out of {}".format(diff,res,fal,unr,cou)























