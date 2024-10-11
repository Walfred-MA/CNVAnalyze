#!/usr/bin/env python3


import os
import collections as cl
import numpy as np
from scipy import stats
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency
import pandas as pd

def find_common(list1, list2):
	
	common = [x for x in list2 if x in list1]
	uncommon2 = [x for x in list2 if x not in list1]
	
	return common, uncommon2

def HLAcomp(list1, list2):
	
	common, uncommon2 = find_common(list1, list2)
	return len(common), len(uncommon2)

def HLAcomp2(list1, list2):
	
	
	list1 = [":".join(x.split(":")[0:3]) for x in list1]
	list2 = [":".join(x.split(":")[0:3]) for x in list2]
	
	common, uncommon2 = find_common(list1, list2)
	return len(common), len(uncommon2)

def HLAcomp3(list1, list2):
	
	
	list1 = [":".join(x.split(":")[0:1]) for x in list1]
	list2 = [":".join(x.split(":")[0:1]) for x in list2]
	
	common, uncommon2 = find_common(list1, list2)
	return len(common), len(uncommon2)


allele_file = "/Users/walfred/Documents/Marklab/benchfigure3/PangenomeAlleles_typefix.tsv_1KG_summary.txt"
immanno_file = "/Users/walfred/Documents/Marklab/HLA/all_Immout"

imm_data = cl.defaultdict(list)
lr_counts = cl.defaultdict(list)
with open(immanno_file, mode = 'r') as f:
	
	for line in f:
		
		line = line.split()
		imm_data[line[0]] = line[1].split(",")
		lr_counts[line[0].split("_")[2]].extend(line[1].split(","))
		
hprc = set(['HG01358', 'HG005', 'HG01109', 'HG00621', 'HG00673', 'HG01106', 'HG03098', 'NA18906', 'HG02145', 'HG01258', 'HG02486', 'HG03540', 'HG01243', 'HG02148', 'HG00741', 'HG02055', 'HG03453', 'HG002', 'HG00733', 'HG03579', 'HG02723', 'HG02559', 'NA21309', 'HG01891', 'HG01978', 'NA20129', 'HG02109', 'NA19240', 'HG03492', 'HG02257', 'HG01928', 'HG03486', 'HG00438', 'HG02630', 'HG02886', 'HG00735', 'HG03516', 'HG02572', 'HG01175', 'HG02818', 'HG01123', 'HG01952', 'HG01361', 'HG01071', 'HG02717', 'HG02622'])


allelenums = [0,0,0,0,0,0,0,0,0,0,0,0]
groups = cl.defaultdict(lambda :[0,0])
total  = 0
signifi = 0
allhwe = []
fullfolder = '/Users/walfred/Documents/Marklab/HLA/T1Ks/'

allfiles = [fullfolder+x for x in os.listdir(fullfolder)]

gp_counts = []

results = cl.defaultdict(lambda: [0,0,0,0,0,0])

samplecount = 0

for file in allfiles:
	
	gp_counts = []
	
	name = file.split("/")[-1].split("_")[0]
	
	if name not in hprc:
		continue
	
	samplecount+=1
	with open(file, mode = 'r') as f:
		
		for line in f:
			
			gp_counts.append(line.split()[0])
	
	
	alltypes = set([x.split("*")[0] for x in lr_counts[name]])
	
	for atype in alltypes:
	
		gp_counts_type = [x for x in gp_counts if x.startswith(atype)]
		lr_counts_type = [x for x in lr_counts[name] if x.startswith(atype)]
	

		tp0, fn0 = HLAcomp(gp_counts_type, lr_counts_type)
		results[atype][0] += tp0
		results[atype][1] += fn0
		
		tp1, fn1 = HLAcomp2(gp_counts_type, lr_counts_type)
		results[atype][2] += tp1
		results[atype][3] += fn1
	
		tp2, fn2 = HLAcomp3(gp_counts_type, lr_counts_type)
		results[atype][4] += tp2
		results[atype][5] += fn2


tp0_all, fn0_all, tp1_all, fn1_all, tp2_all, fn2_all = 0, 0 ,0 , 0 , 0 ,0
for atype in sorted(list(results.keys())):
	
	if atype.startswith("KIR"):
		continue
	
	result = results[atype]
	
	tp0, fn0, tp1, fn1, tp2, fn2 = results[atype]

	print(atype, fn0/(tp0+fn0),  fn1/(tp1+fn1), fn2/(tp2+fn2) )
	

	tp0_all += tp0
	fn0_all += fn0
	tp1_all += tp1
	fn1_all += fn1
	tp2_all += tp2
	fn2_all += fn2
	
print(tp0_all, fn0_all, tp1_all, fn1_all, tp2_all, fn2_all)	