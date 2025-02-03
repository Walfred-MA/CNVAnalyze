#!/usr/bin/env python3


import os
import collections as cl
import numpy as np
from scipy import stats
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency
import pandas as pd

def find_common(list1, list2):
	
	common = [x for x in list1 if x in list2]
	uncommon2 = [x for x in list1 if x not in list2]
	
	return common, uncommon2

def HLAcomp(list1, list2):
	
	common, uncommon2 = find_common(list1, list2)
	return len(common), len(uncommon2)

def HLAcomp2(list1, list2):
	
	
	list1 = [x.split("*")[0]+x.split("*")[1][:5] for x in list1]
	list2 = [x.split("*")[0]+x.split("*")[1][:5] for x in list2]
	
	common, uncommon2 = find_common(list1, list2)
	return len(common), len(uncommon2)

def HLAcomp3(list1, list2):
	
	
	list1 = [x.split("*")[0]+x.split("*")[1][:3] for x in list1]
	list2 = [x.split("*")[0]+x.split("*")[1][:3] for x in list2]
	
	common, uncommon2 = find_common(list1, list2)
	return len(common), len(uncommon2)


allele_file = "/Users/walfred/Documents/Marklab/benchfigure3/PangenomeAlleles_typefix.tsv_1KG_summary.txt"
immanno_file = "/Users/walfred/Documents/Marklab/HLA/all_Immoutfix.txt"
assemanno_file = "/Users/walfred/Documents/Marklab/HLA/all_Immoutassem.txt"

imm_data = cl.defaultdict(list)
lr_counts = cl.defaultdict(list)
with open(immanno_file, mode = 'r') as f:
	
	for line in f:
		
		if line.startswith("KIR") == False:
			continue
		
		
		line = line.split()
		prefix = line[0].split("_")[0]
		genes=line[1].split(",")
		genes = [x for x in genes if x.startswith(prefix)]
		imm_data[line[0]].extend(genes)
		lr_counts[line[0].split("_")[2]].extend(genes)
		
		
alltypes = set([x.split("*")[0] for name in lr_counts.keys() for x in lr_counts[name]])

lr_counts = cl.defaultdict(list)
with open(assemanno_file, mode = 'r') as f:
	
	for line in f:
		
		line = line.split()
		name = line[0].split('#')[0]
		genes=line[1].split(",")
		lr_counts[name].extend(genes)
		
hprc = set(['HG01358', 'HG005', 'HG01109', 'HG00621', 'HG00673', 'HG01106', 'HG03098', 'NA18906', 'HG02145', 'HG01258', 'HG02486', 'HG03540', 'HG01243', 'HG02148', 'HG00741', 'HG02055', 'HG03453', 'HG002', 'HG00733', 'HG03579', 'HG02723', 'HG02559', 'NA21309', 'HG01891', 'HG01978', 'NA20129', 'HG02109', 'NA19240', 'HG03492', 'HG02257', 'HG01928', 'HG03486', 'HG00438', 'HG02630', 'HG02886', 'HG00735', 'HG03516', 'HG02572', 'HG01175', 'HG02818', 'HG01123', 'HG01952', 'HG01361', 'HG01071', 'HG02717', 'HG02622'])


allelenums = [0,0,0,0,0,0,0,0,0,0,0,0]
groups = cl.defaultdict(lambda :[0,0])
total  = 0
signifi = 0
allhwe = []
fullfolder = '/Users/walfred/Documents/Marklab/HLA/fullset2/'

allfiles = [fullfolder+x for x in os.listdir(fullfolder) if x.split(".")[0] in hprc]

results = cl.defaultdict(lambda: [0,0,0,0,0,0])

gp_counts = []

tp0_all, fn0_all, tp1_all, fn1_all, tp2_all, fn2_all = 0, 0 ,0 , 0 , 0 ,0
for file in allfiles:
	
	gp_counts = []
	
	name = file.split("/")[-1].split(".")[0]
	
	with open(file, mode = 'r') as f:
		for line in f:
			
			if line.startswith("result") == False :
				continue
			
			line = line.strip().split()
			
			
			if len(line) < 2 or line[1].startswith("KIR") == False:
				continue
			
			
			genotypes = line[1].split(",")[:-1]
			
			for allele in genotypes:
				
				gp_counts.extend(imm_data[allele])
				
				
				
	for atype in alltypes:
		
		gp_counts_type = [x for x in gp_counts if x.startswith(atype+"*")]
		lr_counts_type = [x for x in lr_counts[name] if x.startswith(atype+"*")]
		
		tp0, fn0 = HLAcomp(gp_counts_type, lr_counts_type)
		results[atype][0] += tp0
		results[atype][1] += fn0
		
		tp1, fn1 = HLAcomp2(gp_counts_type, lr_counts_type)
		results[atype][2] += tp1
		results[atype][3] += fn1
		
		tp2, fn2 = HLAcomp3(gp_counts_type, lr_counts_type)
		results[atype][4] += tp2
		results[atype][5] += fn2
		
		
		tp0_all += tp0
		fn0_all += fn0
		tp1_all += tp1
		fn1_all += fn1
		tp2_all += tp2
		fn2_all += fn2
		
print(tp0_all, fn0_all, tp1_all, fn1_all, tp2_all, fn2_all)
for atype in sorted(list(results.keys())):
	
	result = results[atype]
	
	tp0, fn0, tp1, fn1, tp2, fn2 = results[atype]
	print(atype, fn0/(tp0+fn0),  fn1/(tp1+fn1), fn2/(tp2+fn2) )
	
	
	
print("LOO \n\n")
LOOfolder = '/Users/walfred/Documents/Marklab/HLA/LOO2/'

allfiles = [LOOfolder+x for x in os.listdir(LOOfolder)if x.split("_")[0] in hprc]

gp_counts = []

results = cl.defaultdict(lambda: [0,0,0,0,0,0])
tp0_all, fn0_all, tp1_all, fn1_all, tp2_all, fn2_all = 0, 0 ,0 , 0 , 0 ,0

for file in allfiles:
	
	gp_counts = []
	
	name = file.split("/")[-1].split(".")[0].split("_")[0]
	
	if name not in hprc:
		continue
	
	
	with open(file, mode = 'r') as f:
		
		for line in f:
			
			if line.startswith("result") == False :
				continue
			
			line = line.strip().split()
			
			if len(line) < 2 or line[1].startswith("KIR") == False:
				continue
			
			
			genotypes = line[1].split(",")[:-1]
			
			for allele in genotypes:
				
				gp_counts.extend(imm_data[allele])
				
	for atype in alltypes:
		
		gp_counts_type = [x for x in gp_counts if x.startswith(atype+"*")]
		lr_counts_type = [x for x in lr_counts[name] if x.startswith(atype+"*")]
		
		tp0, fn0 = HLAcomp(gp_counts_type, lr_counts_type)
		results[atype][0] += tp0
		results[atype][1] += fn0
		
		tp1, fn1 = HLAcomp2(gp_counts_type, lr_counts_type)
		results[atype][2] += tp1
		results[atype][3] += fn1
		
		tp2, fn2 = HLAcomp3(gp_counts_type, lr_counts_type)
		results[atype][4] += tp2
		results[atype][5] += fn2
		
		tp0_all += tp0
		fn0_all += fn0
		tp1_all += tp1
		fn1_all += fn1
		tp2_all += tp2
		fn2_all += fn2
		
print(tp0_all, fn0_all, tp1_all, fn1_all, tp2_all, fn2_all)	

for atype in sorted(list(results.keys())):
	
	result = results[atype]
	
	tp0, fn0, tp1, fn1, tp2, fn2 = results[atype]
	
	print(atype, fn0/(tp0+fn0),  fn1/(tp1+fn1), fn2/(tp2+fn2) )
	
	