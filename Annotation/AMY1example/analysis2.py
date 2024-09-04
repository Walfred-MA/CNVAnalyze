#!/usr/bin/env python3
import collections as cl
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from ete3 import Tree

colors = ['lightcyan', 'lightpink', 'lightgreen',  'lightgrey', 'lightskyblue', 'lightsalmon', 'lightslategray',  'lightyellow', 'lightskyblue',   'lightgray', 'lightcoral', 'lightseagreen', 'lightsteelblue', 'lightcyan', 'lightpink',  'lightgreen', 'lightgrey', 'lightskyblue', 'lightsalmon', 'lightslategray', 'lightyellow', 'lightskyblue', 'lightgray',  'lightcoral', 'lightseagreen', 'lightsteelblue', 'lightcyan', 'lightpink', 'lightgreen', 'lightgrey']


thefile = "./PangenomeAlleles_type_AMY.tsv"

table = pd.read_csv(thefile, header= None,sep = '\t')


svtypes = [[0,1,2,3], [4], [5,6,7,8,9,10],[11], [12,13,14,15], [16],[17,18,19,20],[21,22,23,24], [25,26,27,28,29,30,31],[32],[33],[34],[35,36,37],[38],[39]]

singleton_sv_types =set()
typetosvtype = dict()
for i,thetype in enumerate(svtypes):
	for index in thetype:
		typetosvtype[index] = i+1
	if len(thetype) == 1:
		singleton_sv_types.add(thetype[0])
	
allhaplos = set()
haplos_cnv = cl.defaultdict(lambda : [0,0,0])
haplos = cl.defaultdict(list)
haplos_type = cl.defaultdict(list)

major_types_classifies = cl.defaultdict(list)
singletons_type = set()
singletons_sv = set()
for i,row in enumerate(table.values.tolist()):
	
	
	gene = row[3]
	if type(gene) == type(""):
		gene = [x.split(":")[0] for x in row[3].split(";")]
	else:
		gene = []
	
	AMY1,AMY2A,AMY2B = len([x for x in gene if x.startswith("AMY1")]), gene.count('AMY2A'), gene.count('AMY2B')

	theclass = row[1]
	row = row[-1]
	
	
	for gene in row.split(","):
		
		if "chr" in gene:
			haplo = 'HG38'
		elif "NC_0609" in gene:
			haplo = "CHM13"
		else:
			haplo = "_".join(gene.split('_')[2:4])
		
		if typetosvtype.get(i):
			haplos[haplo].append(typetosvtype[i])
			haplos_type[haplo].append(i)
			
			major_types_classifies[typetosvtype[i]].append(theclass)
		
			haplos_cnv[haplo][0]+= AMY1
			haplos_cnv[haplo][1]+= AMY2A
			haplos_cnv[haplo][2]+= AMY2B
			
			if len(row.split(",")) == 1:
				singletons_type.add(haplo)
				
				if i in singleton_sv_types:
					singletons_sv.add(haplo)

print(haplos_type)
haplos = {haplo:sorted(data) for haplo,data in haplos.items()}
haplos_type = {haplo:sorted(data) for haplo,data in haplos_type.items()}

svtypes = cl.defaultdict(list)
for haplo, data in haplos.items():
	svtypes [tuple(sorted(data))].append(haplo)
	
	allhaplos.add(haplo)

svtypes_2 = cl.defaultdict(list)
for haplo, data in haplos_type.items():
	svtypes_2[tuple(sorted(data))].append(haplo)


svtypes_3 = cl.defaultdict(list)
for haplo, data in haplos_cnv.items():
	svtypes_3[tuple(sorted(data))].append(haplo)

singletons = [y for x,y in svtypes.items() if len(y) == 1]
singletons_2 =  [y for x,y in svtypes_2.items() if len(y) == 1]
singletons_3 =  [y for x,y in svtypes_3.items() if len(y) == 1]


#sv types
print(len(singletons))

#subtypes
print(len(singletons_2))

#the aggregate CN
print(len(singletons_3))
print(len(allhaplos))

cnvs = [x for x,y in haplos_cnv.items() if sum(y) > sum(haplos_cnv['HG38']) ]


#high copy numbers

singletons = [y for x,y in svtypes.items() if len(y) == 1 and y[0] in cnvs]
singletons_2 =  [y for x,y in svtypes_2.items() if len(y) == 1 and y[0] in cnvs]
singletons_3 =  [y for x,y in svtypes_3.items() if len(y) == 1 and y[0] in cnvs]


print(len(singletons))
print(len(singletons_2))
print(len(singletons_3))
print(len(cnvs))



singletons = [y for x,y in svtypes.items() if len(y) == 1 and y[0] in singletons_type]
singletons_2 =  [y for x,y in svtypes_2.items() if len(y) == 1 and y[0] in singletons_type]
singletons_3 =  [y for x,y in svtypes_3.items() if len(y) == 1 and y[0] in singletons_type]

print(len(singletons_type))
print(len(singletons_sv))
print(len(singletons))
print(len(singletons_2))
print(len(singletons_3))