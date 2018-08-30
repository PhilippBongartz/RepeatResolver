#!/usr/bin/env python
# coding: latin1

import random
import math
import sys

# This script creates simulated data to test the repeatresolver.

# Sub 1.4, Ins 11.5, Del 3.4, Match 83.7
NotBase={'a':['c','g','t'],'c':['a','g','t'],'g':['c','a','t'],'t':['c','g','a']}
def PacBioError(seq):
	read=''
	for x in range(len(seq)):
		rand=random.random()
		if rand<0.837+0.115:
			read+=seq[x]
		elif rand<0.837+0.115+0.014:
			read+=NotBase[seq[x]][int(random.random()*3)]
		elif rand<0.837+0.115+0.014+0.034:
			pass
		# Insertion separately otherwise there will always be just one insertion:
		rand=random.random()
		while rand<0.103139:    # geometric formula: 1.0/(1.0-0.103139) - 1.0 = 0.115
			read+='acgt'[int(random.random()*4)]
			rand=random.random()
	return read

def DistributedVarCopies(seq,copynumber,difference):
	SNPnumber=int(len(seq)*difference*3)
	Copies=[seq for t in range(copynumber)]
	Positions=[10+int(random.random()*(len(seq)-20)) for t in range(SNPnumber)]
	Positions.sort()
	for t in range(SNPnumber):
		position=Positions[-t-1]
		random.shuffle(Copies)
		rand=int(random.random()*len(Copies)) # Every SNP on another subset
		errortype=random.random()
		if errortype<=1.0/3.0: # Sub
			for x in range(rand):
				Copies[x]=Copies[x][:position]+NotBase[Copies[x][position]][rand%3]+Copies[x][position+1:]
		if 1.0/3.0<errortype<=2.0/3.0: # Del
			for x in range(rand):
				Copies[x]=Copies[x][:position]+Copies[x][position+1:]			
		if 2.0/3.0<errortype<=3.0/3.0: # Ins
			base='acgt'[int(4*random.random())]
			for x in range(rand):
				Copies[x]=Copies[x][:position]+base+Copies[x][position:]
	return Copies	


def GradientCopies(seq,copynumber,difference):
	Copies=[seq]
	SNPnumber=int(len(seq)*difference)
	for c in range(copynumber-1):
		copy=Copies[-1]
		for t in range(SNPnumber):
			position=int(random.random()*len(copy))
			errortype=random.random()
			rand=int(random.random()*3)
			if errortype<=1.0/3.0: # Sub
				copy=copy[:position]+NotBase[copy[position]][rand]+copy[position+1:]
			if 1.0/3.0<errortype<=2.0/3.0: # Del
				copy=copy[:position]+copy[position+1:]			
			if 2.0/3.0<errortype<=3.0/3.0: # Ins
				base='acgt'[int(4*random.random())]
				copy=copy[:position]+base+copy[position:]
		Copies.append(copy)
	return Copies


def EquiDistantRepeatCopies(seq,copynumber,difference):
	difference/=2.0 # If each copy gets .5 diff they differ by 1.
	SNPnumber=int(difference*len(seq))
	Copies=[]
	for c in range(copynumber):
		copy=seq
		for t in range(SNPnumber):
			position=int(random.random()*len(seq))
			errortype=random.random()
			rand=int(random.random()*3)
			if errortype<=1.0/3.0: # Sub
				copy=copy[:position]+NotBase[copy[position]][rand]+copy[position+1:]
			if 1.0/3.0<errortype<=2.0/3.0: # Del
				copy=copy[:position]+copy[position+1:]			
			if 2.0/3.0<errortype<=3.0/3.0: # Ins
				base='acgt'[int(4*random.random())]
				copy=copy[:position]+base+copy[position:]
		Copies.append(copy)
	return Copies


def TreeCopies(seq,copynumber,difference):
	difference/=2.0 
	SNPnumber=int(difference*len(seq))
	Copies=[[seq],[]]
	for t in range(int(math.log(copynumber,2))+1):
		for oldcopy in Copies[t%2]:
			for c in range(2):
				copy=oldcopy
				for tt in range(SNPnumber):
					position=int(random.random()*(len(oldcopy)-SNPnumber))
					errortype=random.random()
					rand=int(random.random()*3)
					if errortype<=1.0/3.0: # Sub
						copy=copy[:position]+NotBase[copy[position]][rand]+copy[position+1:]
					if 1.0/3.0<errortype<=2.0/3.0: # Del
						copy=copy[:position]+copy[position+1:]			
					if 2.0/3.0<errortype<=3.0/3.0: # Ins
						base='acgt'[int(4*random.random())]
						copy=copy[:position]+base+copy[position:]

				Copies[(t+1)%2].append(copy)
		Copies[t%2]=[]	
	return Copies[(t+1)%2][:copynumber]



def RandomSequence(length):
	seq=''
	for t in range(length):
		seq+='acgt'[int(4*random.random())]
	return seq

# Length distribution taken from the Drosophila Reads:
Lengthshisto=[0, 323, 427, 411, 355, 353, 358, 321, 293, 321, 281, 275, 241, 239, 226, 
	185, 177, 162, 126, 117, 126, 108, 88, 83, 61, 52, 51, 29, 16, 7, 3, 1, 1, 0, 0, 0, 0, 0, 0, 0]

# Sampling reads according to the length distribution of the histone reads from a given complex sequence
def ReadSampling(coverage, Lengthshisto, genome):
	Lengthsprob=[float(Lengthshisto[t])/float(sum(Lengthshisto)) for t in range(len(Lengthshisto))]
	Lengths=[]
	CovLengths=[]
	Starts=[]
	current_coverage=0
	while current_coverage<coverage:
		rand=random.random()
		length=-1
		prob=0
		while prob<rand:
			length+=1
			prob+=Lengthsprob[length]

		length*=1000
		length+=int(random.random()*1000)
		Lengths.append(length)
		start=int(random.random()*(len(genome)-length))
		Starts.append(start)

		covlength=min(len(genome)-10000,start+length)-max(start,10000)  # the length of the repetitive part
		CovLengths.append(covlength)    
		current_coverage=float(sum(CovLengths))/float(len(genome)-20000)

	Reads=[]
	for l in range(len(Lengths)):
		start=Starts[l]
		length=Lengths[l]
		Reads.append(PacBioError(genome[start:start+length]))

	return Reads,Starts

# parameters

coverage=40
copynumber=100
difference=0.01
repeatlength=30000
Type='Tree'

# Reading in Parameters: 
if len(sys.argv)>1:
	for t in range(1,len(sys.argv)-1):
		if sys.argv[t]=='-c':
			coverage=int(sys.argv[t+1])
		if sys.argv[t]=='-n':
			copynumber=int(sys.argv[t+1])
		if sys.argv[t]=='-d':
			difference=float(sys.argv[t+1])/100
		if sys.argv[t]=='-l':
			repeatlength=int(sys.argv[t+1])
		if sys.argv[t]=='-t':
			Type=sys.argv[t+1]

# Check:
if Type not in ['Distributed','EquiDistant','Tree']:
	print "The type {} is not correct: Distributed, EquiDistant or Tree.".format(Type)
	sys.exit()

# removing '.' and trailing zeroes to not mess up the file name too much. 
percstring=str(difference*100.0)
i=len(percstring)-1
while percstring[i]=='0':
	percstring=percstring[:i]
	i-=1
if percstring[i]=='.':
	percstring=percstring[:i]

datasetname=Type+'_'+percstring.replace('.','')+'perc_'+str(repeatlength)+'kb'
print "The data set name:",
print datasetname



# Die Repeat Sequence

seq=RandomSequence(repeatlength)

# The Copies

if Type=='Tree':
	Copies=TreeCopies(seq,copynumber,difference)

if Type=='Distributed':
	Copies=DistributedVarCopies(seq,copynumber,difference)

if Type=='EquiDistant':
	Copies=EquiDistantRepeatCopies(seq,copynumber,difference)

#Copies=GradientCopies(seq,copynumber,difference)

# Adding flanking sequences:
for c in range(len(Copies)):
	left=RandomSequence(10000)
	right=RandomSequence(10000)
	Copies[c]=left+Copies[c]+right

# ReadSampling from all Copies
AllReads=[]
AllInfo=[]
for c in range(len(Copies)):
	Reads,Starts=ReadSampling(coverage,Lengthshisto, Copies[c])
	AllReads+=Reads 
	AllInfo+=[(s,c) for s in Starts]

print "Reads: {}".format(len(AllReads))

# Output Reads and Start+Copy Information


# The simulated reads sampled with pacbio error from the simulated complex
f=open(datasetname+'.fasta','w')
for read in AllReads:
	f.write('>\n')
	for t in range(0,len(read),100):
		f.write(read[t:t+100]+'\n')
f.close()

# The correct placements of the reads as ground truth
f=open(datasetname+'_ReadPlacements','w')
for start,copy in AllInfo:
	f.write(str(start)+'\n')
f.close()

f=open(datasetname+'_ReadCopynumbers','w')
for start,copy in AllInfo:
	f.write(str(copy)+'\n')
f.close()

f=open(datasetname+'_Template.fasta','w')
f.write('>\n')
f.write(seq+'\n')
f.close()




