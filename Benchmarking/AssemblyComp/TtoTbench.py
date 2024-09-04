#!/usr/bin/env python3


import collections as cl

typefiles = "/Users/walfred/Documents/Marklab/benchfigure3/PangenomeAlleles_typefix.tsv"

CNgenotype =  "/Users/walfred/Documents/Marklab/benchfigure3/CN00001_out_leave2.txt_allelefix.txt"
assemname = "CN00001"

def category(x):
	
	if x <= 1:
		return "=1"
	elif x <=2:
		return "=2"
	elif x<=5:
		return "<=5"
	elif x<10:
		return "<10"
	elif x<20:
		return "<20"
	else:
		return ">=20"
	
allele_types = cl.defaultdict(list)
geneto_types = dict()

groupcount_lr = cl.defaultdict(int)
genotypecount_lr =  cl.defaultdict(int)

singletons = set()

filtername = set([])
chm13_counts = cl.defaultdict(int)

typesmember = cl.defaultdict(list)

decoys = set()

with open(typefiles, mode = 'r') as f:
	
	for line in f:
		if len(line) == 0:
			continue
		
		line =  line.strip().split('\t')
		allelename =  line[0]
		groupname = "_".join(allelename.split('_')[:2])
		
		#no decoy or  "NA" in line[3]
		if line[5] in ["Intron","Decoy"] :
			decoys.add(allelename)
			
			
			
with open(typefiles, mode = 'r') as f:
	
	for line in f:
		if len(line) == 0:
			continue
		"""
		if line.split()[3] == '0' or int(line.split()[5]) <1000:
			filtername.add(line.split()[0])
			continue
		"""
		
		line =  line.strip().split('\t')
		
		allelename =  line[0]
		
		groupname = "_".join(allelename.split('_')[:2])
		
		kmernum = int(line[6])
		
		genes = [x.split("@")[0] for x in line[-1].split(",") if ":" not in x]
		
		allele_types[allelename] = genes
		
		for gene in genes:
			geneto_types[gene] = line[0]
			
			
		chm13_counts[groupname] += len([x for x in genes if "NC_" in x])
		
		lrgenes = [x for x in genes if assemname in x]
		lrcount = len(lrgenes)
		
		#no decoy 
		if allelename in decoys:
			continue
		
		groupcount_lr[groupname] += lrcount
		genotypecount_lr[allelename] = lrcount
		
				
		if lrcount == len(genes):
			
			singletons.add(allelename)
			
			
groupcount = cl.defaultdict(int)
genotypecount = cl.defaultdict(int)
with open(CNgenotype, mode = 'r') as f:
	for line in f:
		
		if len(line.strip()) == 0 or ":" in line or ">" in line:
			continue
		line = line.strip()
		
		allelename = line.split()[0]
		
		if allelename in decoys:
			continue
		
		
		groupname = "_".join(allelename.split('_')[:2])
		
		coefs = line.split()[-1].split(",")
		coefs = sum([float(x) for x in coefs])
		
		groupcount[groupname] += int(coefs+0.5)
		
		genotypecount[allelename] += int(coefs+0.5)
		
#filtered decoys 
for typename, count in groupcount_lr.items():
	
	if typename not in groupcount and count > 0 :
		
		groupcount[typename] = 0
		
for typename, count in genotypecount_lr.items():
	
	if typename not in genotypecount and count > 0 :
		
		genotypecount[typename] = 0

for typename, count in groupcount.items():
	
	if typename not in groupcount_lr and count > 0 :
		
		groupcount_lr[typename] = 0
		
for typename, count in genotypecount.items():
	
	if typename not in genotypecount_lr and count > 0 :
		
		genotypecount_lr[typename] = 0
		
		
genosum1 = 0
assemsum1 = 0
total1 = 0
error1 = 0
aggrecomp = []

bins = cl.defaultdict(lambda: [0,0,0,0,0])
for typename, count in groupcount.items():
	
	if count < 0:
		continue
	
	if typename in groupcount_lr:
		
		cate = category(chm13_counts[typename])
		
		genosum1 += count
		assemsum1 += groupcount_lr[typename]
		
		total1 += count
		error1 += abs(count-groupcount_lr[typename])
		
		cate = category(chm13_counts[typename])
		bins[cate][0] += min(count, groupcount_lr[typename])
		bins[cate][1] += max(0,count-groupcount_lr[typename])
		bins[cate][2] += max(0,groupcount_lr[typename]-count)
		
		aggrecomp.append([chm13_counts[typename],count, groupcount_lr[typename]])
		
		print("aggregate:",typename,chm13_counts[typename],count, groupcount_lr[typename])
		
		
typeerror = cl.defaultdict(lambda: [0,0])
genosum2 = 0
assemsum2 = 0
error2 = 0
typecomp = []

errors = []
for name, count in genotypecount.items():
	
	
	if name in genotypecount_lr:
		
		groupname = "_".join(name.split('_')[:2])
		genosum2 += min(count,genotypecount_lr[name])
		assemsum2 += genotypecount_lr[name]
		error2 += abs(count-genotypecount_lr[name])
		
		cate = category(chm13_counts[groupname])
		bins[cate][3] += abs(count-genotypecount_lr[name])
		
		if abs(count-genotypecount_lr[name]):
			errors.extend(typesmember[name])
			print(name, count, genotypecount_lr[name])
			
		typeerror[groupname][0] += min(count, genotypecount_lr[name] )
		typeerror[groupname][1] += abs(count-genotypecount_lr[name])
		
		if name in singletons:
			
			bins[cate][4] += 1
			
		typecomp.append([chm13_counts[groupname], count, genotypecount_lr[name]])
		
		
print(genosum1 , assemsum1, total1, error1)
print(genosum2 , assemsum2, error2)

binnames = sorted(list(bins.keys()))

for binname, bin in bins.items():
	
	printline = list(map(str,[binname]+bin))
	
	print(",".join(printline ))
	
	
keysort = sorted(list(typeerror.keys()), key = lambda x: typeerror[x][1]/max(1,typeerror[x][0]), reverse =1)

