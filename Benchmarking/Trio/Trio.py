#!/usr/bin/env python3



import collections as cl
import numpy as np
from scipy import stats
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency
import pandas as pd

thefile = "./integrated_call_samples_v3.20200731.ALL.ped"
	
table = pd.read_csv(thefile,sep = "\t")

parents = dict()

unrelates = []
for index, row in table.iterrows():
	
	if row['Paternal ID'] != '0' and row['Maternal ID'] != '0':
	
		parents[row['Individual ID']] = [row['Paternal ID'], row['Maternal ID']]
	
	
probands = list(parents.keys())


allele_file = "./PangenomeAlleles_typefix.tsv_1KG_summary.txt"

trio_out = allele_file+"_trio.out"

bench_posi = [0,0,0]
bench_nega = [0,0,0]
bench_mix = [0,0,0]
bench_dup = [0,0,0]

bench_fp = cl.defaultdict(lambda: [0,0])
bench_fn = cl.defaultdict(lambda: [0,0])

groups = cl.defaultdict(lambda: [0,0])

samplecorr = 0

hprc = set(['HG002', 'HG00438', 'HG005', 'HG00621', 'HG00673', 'HG00733', 'HG00735', 'HG00741', 'HG01071', 'HG01106', 'HG01109', 'HG01123', 'HG01175', 'HG01243', 'HG01258', 'HG01358', 'HG01361', 'HG01891', 'HG01928', 'HG01952', 'HG01978', 'HG02055', 'HG02080', 'HG02109', 'HG02145', 'HG02148', 'HG02257', 'HG02486', 'HG02559', 'HG02572', 'HG02622', 'HG02630', 'HG02717', 'HG02723', 'HG02818', 'HG02886', 'HG03098', 'HG03453', 'HG03486', 'HG03492', 'HG03516', 'HG03540', 'HG03579', 'NA18906', 'NA19240', 'NA20129', 'NA21309'])

allelenums = [0,0,0,0,0,0,0,0,0,0,0,0]
groups = cl.defaultdict(lambda :[0,0])
total  = 0
signifi = 0
allhwe = []
with open(allele_file, mode = 'r') as f:
	
	sample_names = [x.split(".")[0] for x in f.readline().split()[1:]]
			
	probands = [x for x in probands if x in sample_names and parents[x][0] in sample_names and parents[x][1] in sample_names ]		
	
	lineindex = 0
	for line in f:
		
		lineindex+=1
		
				
		line = line.strip().split('\t')
		
		name = line[0]
		allhaplos = ["_".join(x.split("_")[2:4]) for x in line[9].split(",") if "_h" in x]
		ifcollapse = len(allhaplos) - len(set(allhaplos))
		ifsex = int ("chrX" in line[4] or "chrY" in line[4])
		ifintronordecoy = int(line[5] != 'Exon')
		groupname = name.split('_')[0]		
		line = line[10:]
		
		allcounts = list(map(float, line))
		
		#len([x for x in allcounts if x >2]) > len(allcounts )
		#name, ifcollapse, ifsex, ifintron ,samplesize
		
		if ifintronordecoy :
			
			continue
			
		sample_floats = {name:count for name,count in zip(sample_names,allcounts)}
		sample_counts = {name:int(count+0.5) for name,count in zip(sample_names,allcounts)}
		
		allcounts = list(sample_counts.values())
		
		
		max_count = max(allcounts)
		
		group = "_".join(name.split("_")[:3])
		name = "_".join(name.split('_')[:2])
		for proband in probands:
			
			proband_count = sample_counts[proband]
			
			parents_name = parents[proband]
			
			parents_count = [sample_counts[parents_name[0]], sample_counts[parents_name[1]]]
			
			known_posi = len([x for x in parents_count if x > 1])
			known_nega =  len([x for x in parents_count if x == 0])
			
			proband_count0 = max(0, 2 - proband_count)
			
			error_posi = max(known_posi  - proband_count, 0 )
			error_nega = max(known_nega - proband_count0 , 0 )
			
			if known_posi :
				bench_fp[name][0] += known_posi
				bench_fp[name][1] += error_posi
				
			if known_nega:
				bench_fn[name][0] += known_nega
				bench_fn[name][1] += error_nega
					
			if proband_count == 0:
				
				bench_posi[0] += 1
				#groups[name][0] += 1
				error = len([x for x in parents_count  if x > 1 ])
				groups[name][0] += error
				if error > 0:
					
					groups[name][1] += error
					bench_posi [error] += error
			
			elif proband_count > 2 :
				
				error =  max(proband_count - sum(parents_count), 0)
				
				bench_dup[0] += proband_count 
				groups[name][0] += 2 * proband_count + error
				if error > 0:
					groups[name][1] += error
					bench_dup [1] += error

			elif proband_count == 2:
				
				error =  max(proband_count - sum([1 for x in parents_count if x > 0]), 0)
				#error =  max(proband_count - sum(parents_count), 0)
				
				bench_nega[0] += proband_count
				groups[name][0] += 2* proband_count + error
				if error > 0:
					groups[name][1] += error
					bench_nega [1] += error
			
			elif proband_count == 1 :
				
				error = 0
				if len([x for x in parents_count if x > 1]) == 2 or len([x for x in parents_count if x == 0]) == 2:
					error =  1
				
				bench_mix[0] += proband_count 
				groups[name][0] += 2 *proband_count + error
				if error > 0:
					groups[name][1] += error
					bench_mix [error] += error
			
			if proband_count > 0:
				pass
				#print(group, proband, parents_name, [proband_count]+parents_count, sample_floats[proband], [sample_floats[x] for x in parents_name])
				
			
			if error > 0 and proband in hprc:
				#print(group, proband, parents_name, [proband_count]+parents_count, sample_floats[proband], [sample_floats[x] for x in parents_name])
				
				pass
				if proband in hprc or parents_name[0] in hprc or parents_name[1] in hprc:
					#print("HPRC")
					pass
				

print(bench_posi)
print(bench_nega)
print(bench_mix)
print(bench_dup)


groups_sort = sorted([(value[1]/max(1,value[0]),value[1],value[0],name) for name, value in groups.items()],reverse = 1)

with open(trio_out,mode = 'w') as f:

	for errate,error, total, name in groups_sort:
		
		f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(errate, error, total,bench_fp[name][0],bench_fp[name][1],bench_fn[name][0],bench_fn[name][1],name))