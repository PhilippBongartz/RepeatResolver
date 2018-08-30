#!/usr/bin/env python
# coding: latin1

import random
import math
import os
import sys
import time




# Here we read a MSA and output several options how to split it into sections that are separately resolved. 

if len(sys.argv)<2:
	print "Missing msa path."
	sys.exit()

path=sys.argv[1]
if os.path.exists(path):
	f=open(path,'r')
	MA=[]
	for line in f.readlines():
		MA.append(line[:-1])
	f.close()

else:
	print "Wrong msa path."
	sys.exit()

# Reading in Parameters: 
coverage=35
parts=6
if len(sys.argv)>1:
	for t in range(1,len(sys.argv)-1):
		if sys.argv[t]=='-c':
			coverage=int(sys.argv[t+1])
		if sys.argv[t]=='-p':
			parts=int(sys.argv[t+1])

Coverages=[sum([1 for z in range(len(MA)) if MA[z][c]!=' ']) for c in range(0,len(MA[0]),100)]

average=sum(Coverages)/len(Coverages)
maximus=max(Coverages)

start=0
while Coverages[start]<coverage:
	start+=1
start*=100

ende=len(Coverages)-1
while Coverages[ende]<coverage:
	ende-=1
ende*=100

print "Subdivision of the MSA into section:"
print start,
for p in range(parts):
	print start+(p+1)*(ende-start)/parts,
print 

sys.exit()
