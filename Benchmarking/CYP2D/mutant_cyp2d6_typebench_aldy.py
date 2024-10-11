#!/usr/bin/env python3

import os
import re
import collections as cl


typefile = "/Users/walfred/Documents/Marklab/HLA/CYPD26.group_table.tsv"

samplecount = cl.defaultdict(set)
haplotypecount  = cl.defaultdict(set)
nametostartype = cl.defaultdict(list)
nametoshaplotype = dict()
with open(typefile, mode = 'r') as f:
	for line in f:
		line = line.split()
		if line[-1] == "None":	
			nametostartype[line[0].split("_h")[0]].append("*0")
		else :
			nametostartype[line[0].split("_h")[0]].append(line[-1])
		
		nametoshaplotype[line[0]] = (line[1],line[-1])
		samplecount[line[-1]].add(line[0].split("_")[0])
		haplotypecount[(line[1],line[-1])].add(line[0].split("_")[0])
		
allfiles = ["/Users/walfred/Documents/Marklab/figure4/results/"+x for x in os.listdir("/Users/walfred/Documents/Marklab/figure4/results")]


TP = 0
FN = 0
FP = 0
nametotypes = cl.defaultdict(list)
for afile in allfiles:
	
	samplename = afile.split("/")[-1].split(".")[0]
	
		
	with open(afile, mode = 'r') as f:
		
		for line in f:
			if line.startswith("#"):
				continue
			line = line.split()
			
			types = re.findall(r'\*\d+',line[3])
			
			nametotypes[samplename ] = types
			break
	TP += len([x for x in nametostartype[samplename] if x in types])
	FN += len([x for x in nametostartype[samplename] if x not in types])
	FP += len([x for x in types if x not in nametostartype[samplename]])
	
		
print(TP,FN,FP-4)
print((FP-4+FN)/(2*TP+FP-4+FN))