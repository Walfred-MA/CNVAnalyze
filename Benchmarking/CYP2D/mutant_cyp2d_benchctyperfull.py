#!/usr/bin/env python3

#/Users/walfred/anaconda3/bin/python3 $filename -i /Users/walfred/Documents/Marklab/figure4/CYP2D_group1_summary.txt -t /Users/walfred/Documents/Marklab/figure4/CYP2D_group1_CYP2D6OOOCYP2D7OOOCYP2D8P.fasta_annotate.fa_types.txt -g /Users/walfred/Documents/Marklab/figure4/CYP2D.gff3 -o ./CYP2D6_bench.png

import re
import collections as cl
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.patheffects as pe
import os 
import argparse
import colorsys


def getcigarsize(text):
	
	allcigars = sum([int(x[:-1]) for x in re.findall(r'\d+[MIX=]', text)])
	
	return allcigars

def maptoseg(allexons, segments):
	
	l = len( segments)
	
	allexons = [x if i %2 == 0 else x -1 for i,x in enumerate(allexons)]
	
	
	allposi =  segments + allexons 
	
	
	allposi_sortindex = sorted(range(len(allposi)), key = lambda i: allposi[i])
	
	overlaps = [[] for x in segments] 
	overlap = overlaps[0]
	seg_start = 0
	for rank, index in enumerate( allposi_sortindex):
		
		posi = allposi[index]
		if index < l :
			overlap = overlaps[index]
			seg_start = allposi[index]
		else:
			overlap.append(index)
			
	overlaps = overlaps[:-1]
	for i,overlap in enumerate(overlaps):
		
		overlap_ = [allposi[x]- allposi[i] if x % 2 == 0 else allposi[x]- allposi[i]+1 for x in overlap]
		
		if len(overlap) % 2 == 1:
			if (overlap[0] - l) % 2:
				overlaps[i] = [0] + overlap_
			else:
				overlaps[i] = overlap_ + [allposi[i+1] - allposi[i]]
		else:
			overlaps[i] = overlap_
			
	return overlaps

def overlapexons(allexons):
	
	if len(allexons) and type(allexons[0]) == type([]):
		l = len(allexons)
		allposi = sum(allexons, [])
	else:
		l = len(allexons)//2
		allposi = allexons
		
	allposi_sortindex = sorted(range(len(allposi)), key = lambda i: allposi[i])
	
	allgroups = []
	current_group = []
	last_endpoint = 0
	overlap_count = 0
	
	for rank, index in enumerate( allposi_sortindex):
		
		posi = allposi[index]
		if index %2 :
			overlap_count -= 1
			
			if overlap_count == 0:
				
				last_endpoint = allposi[index]
				current_group.append(posi)
		else:
			
			overlap_count += 1
			
			if overlap_count == 1 :
				
				allgroups.append(current_group)
				current_group = [posi]
				
	allgroups.append(current_group)
	
	return allgroups[1:]

def cigar_findrposi(cigar, find_qposes):
	
		l = len(find_qposes)
		if l == 0:
			return []
	
		curr_rposi = 0
		curr_qposi = 0
		new_rposi = 0
		new_qposi = 0
	
		curr_size = 0
	
		find_qpos_index = 0
		find_qpos = find_qposes[0]
		find_rposes = []
		find_qposes = list(find_qposes)+[-1]
	
	
		for char in cigar:
			
				if char <= '9':
					
						curr_size *= 10
						curr_size += ord(char) - ord('0')
						continue
			
				if char == '=' or char ==  'M':
					
						new_rposi = curr_rposi + curr_size
						new_qposi = curr_qposi + curr_size
					
				elif char == 'I':
					
						new_qposi = curr_qposi + curr_size
					
				elif char == 'D' or char == 'H':
					
						new_rposi = curr_rposi + curr_size
					
				elif char == 'X':
					
					
						new_rposi = curr_rposi + curr_size
						new_qposi = curr_qposi + curr_size
					
					
				elif char == 'S':
					
						if curr_qposi == 0:
								curr_qposi = curr_size
							
						curr_size = 0
						continue
			
				while find_qpos_index<l and find_qpos < new_qposi:
					
						if char == '=' or char ==  'M' or char =='X' or char == 'N':
								find_rposes.append(curr_rposi+(find_qpos-curr_qposi))
						else:
								find_rposes.append(curr_rposi)
							
						find_qpos_index += 1
						find_qpos = find_qposes[find_qpos_index]
					
					
				curr_size = 0
				curr_rposi = new_rposi
				curr_qposi = new_qposi
			
		while find_qpos_index<l:
			
				find_qpos_index += 1
				find_rposes.append(curr_rposi)
			
		return find_rposes

def adjust_lightness(color, amount=0.5):
	try:
		c = mcolors.cnames[color]
	except:
		c = color
	c = colorsys.rgb_to_hls(*mcolors.to_rgb(c))
	return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


def segfiltration(cigar, usednames = ""):
	
	cigars = re.findall('[<>][0-9_A-Za-z:]+',cigar)

	newcigars = []
	for cigarstr in cigars:
		
			thename,thecigar = cigarstr.split(":")
			if type(usednames) != type("") and thename[1:] not in usednames:
				
					continue
		
			thesize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=]', thecigar)]+[0])
		
			if thesize >0 and thesize < 50:
				
					fullsize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', thecigar)]+[0])
				
					thecigar = "{}H".format(fullsize )
				
			newcigars.append(thename+":"+thecigar)
		
	newcigars= "".join(newcigars)

	return newcigars

def getgenelocation(table, gff3file, used_cigars):
		
	segment_exonalign = cl.defaultdict(lambda : cl.defaultdict(list))
	gffrecord = cl.defaultdict(lambda : cl.defaultdict(list))
	if len(gff3file):
		
		with open(gff3file, mode = 'r') as f:
			for line in f:
				line = line.strip().split()
				
				if line[2] == 'exon':
					genename = [x for x in line[-1].split(";") if x.startswith("gene_name=")][0][10:]
					
					gffrecord[line[0]][genename].append([int(line[3]), int(line[4])])
						
		gffrecords = {thechr:{gene:overlapexons(coordi) for gene, coordi in record.items()} for thechr, record in gffrecord.items()}
		
		for row in table.values.tolist():
			
			locus = row[7]
			refchr = locus.split(":")[0]
			if refchr not in gffrecords:
				continue
			
			strd = locus[-1]
			
			thepath = row[-3]
			
			cigar = row[-2]
			segment_cigars = re.findall('[<>][0-9_A-Za-z:]+',cigar)
			segment_orders = re.findall('[<>][0-9_A-Za-z:]+',thepath)
			
			segment_cigars = {x.split(":")[0][1:]:x for x in segment_cigars if x.split(":")[0][1:].split("_")[0] in used_cigars}
			segment_orders = [x for x in segment_orders if x[1:].split("_")[0] in used_cigars]
			
			segment_cigars = [segment_cigars[x[1:]] for x in segment_orders]
			
			segment_alignsizes = [(x.split(":")[0], getcigarsize(x.split(":")[1])) for x in segment_cigars ]
			
			qcoordi = 0
			segment_qranges = [0]
			segment_names = []
			for segment in segment_alignsizes:
				segment_names.append(segment[0])
				segment_qranges.append(qcoordi + segment[1])
				qcoordi += segment[1]
			segment_qranges.append(qcoordi)
			
			start, end = locus.split(":")[1][:-1].split('-')
			start, end = int(start), int(end)
			
			for gene, gffrecord in gffrecords[refchr].items():
				
				exons = [x for x in gffrecord if max(x[1],end) - min(x[0],start) - (end - start) - (x[1]-x[0]) < 0 ] 
				
				if len(exons) == 0:
					continue
				
				if strd == '+' :
					exons = sum([ [ max(0, x[0] - start), x[1] - start ] for x in exons],[])
				else:
					exons = sum([ [ max(0, end - 1 - x[1]), end - x[0] ] for x in exons[::-1]],[])
					
				exons_onsegments = maptoseg(exons, segment_qranges)
				
				
				for cigar, exons_onsegment, size in zip(segment_cigars, exons_onsegments, segment_alignsizes):
					
					if len(exons_onsegment) == 0:
						continue
					
					name, cigar = cigar.split(":")
					
					if name[0] == '<':
						exons_onsegment = [size[1]-x for x in exons_onsegment]
						
					exons_onsegment_rposi = cigar_findrposi(cigar, exons_onsegment)
					
					segment_exonalign[gene][name[1:]].extend(exons_onsegment_rposi)
					
		for gene, segment_exonalign_ in segment_exonalign.items():	
			
			for name, exonalign in segment_exonalign_.items():
				
				exonalign = overlapexons(exonalign)
				
				exonalign = [x for x in exonalign if x[1]-x[0] > 50]
				
				segment_exonalign[gene][name] = exonalign
				
	return segment_exonalign

	
def plotexons(fig, ax, lelement, inputfile, gff3file, used_cigars):
	
	table = pd.read_csv(inputfile, sep = '\t', header = None)
	
	row0 = table.values.tolist()[0]
	segment_cigars = re.findall('[<>][0-9_A-Za-z:]+',row0[-2])
	
	segment_order = []
	for segment_cigar in segment_cigars:
			segment_order.append(segment_cigar.split(':')[0][1:])
		
	exonlocations = getgenelocation(table, gff3file, used_cigars)
		
	segment_order = [x for x in segment_order if x in used_cigars]
	exonlocations_order = {gene: [exonlocation[segment.split('_')[0]] for segment in segment_order]  for gene,exonlocation in exonlocations.items() }
	
	max_dupnum = max([int(x.split('_')[-1]) if '_' in x else 1 for x in segment_order])
	
	c = 'red'
	l = 1
	L= 2
	chromy = lelement + 5
	fullsize = sum(list(used_cigars.values()))
	for dupnum in [0,1]:
	
		for gene, exonlocation in exonlocations_order.items():
			
			if len(sum(exonlocation, [])) == 0:
				continue
			
			leftend = 0
			exons_posi = []
			for name, pathexon in zip(segment_order, exonlocation):
				
				dupnum_ = 1 if "_" in name else 0
				
				size = used_cigars[name.split('_')[0]]
				rightend = leftend + size
				
				if dupnum_ == dupnum :
					for exon in pathexon:
						exons_posi.append([ (leftend + exon[0]) / fullsize, (leftend + exon[1]) / fullsize])
				leftend = rightend
			
			exons_posi_ = sum(exons_posi, [])
			
			if len(exons_posi_) == 0:
				continue
			
			print(dupnum, gene)
			
			exons_posi_min, exons_posi_max = min(exons_posi_), max(exons_posi_)
			ax.plot([exons_posi_min, exons_posi_max], [chromy, chromy], color=c, linewidth=l, linestyle='-')
			for exon in exons_posi:
				
				ax.plot(exon, [chromy, chromy], color=c, linewidth=L, linestyle='-')
				
			chromy += 1
	
	return 

"""
def mergesegments(cigar):

		cigars = re.findall(r'[><0-9_:]*\d+[HMX=ID]', cigarstr)

		current_type = ''
		current_size = 0
		name = ""
		newname = ""
		names = []
		newcigars = []
		for cigar in cigars:

				thetype = cigar[-1]
				size = int(cigar.split(":")[-1][:-1])

				if ":" in cigar:
						newname = cigar.split(":")[0].split("_")[0][1:]

				if thetype != current_type or newname != name:
						newcigars.append("{}{}".format(current_size, current_type))
						names.append(name)
						current_size = 0

				current_size += size
				current_type = thetype

		newcigars.append("{}{}".format(current_size, thetype))
		names.append(name)

		newcigars = newcigars[1:]

"""

def getsegments(cigarstr,highlight,exons):
	
	cigars = re.findall(r'[><0-9_:]*\d+[HMX=ID][ATCGNatcgn]*', cigarstr)
	variants = []
	segments = []
	highlight_points = {}
	exon_points = {}
	
	highlight_coordinate = -1
	highlight_type = cl.defaultdict(int)
	
	oldname = ""
	variant = []
	
	localpath_offsite = 0
	coordinate = 0
	laststart = 0
	findposition = -1
	
	for cigar in  cigars:
		
		if ":" in cigar:
			
			localpath_offsite = coordinate
			
			name,cigar = cigar.split(":")
						
			name = name[1:]
			
			if name in highlight:
				
				highlight_points[name] = [[x + coordinate for x in seg] for seg in highlight[name] ]
				
			if name in exons:
				
				exon_points[name] =  [[x + coordinate for x in seg] for seg in exons[name] ]
				
		lcigar = len(re.findall(r'\d+', cigar)[0])
		cigar,seq = cigar[:(lcigar+1)],cigar[(lcigar+1):]
		
		ifmask = seq.islower() if len(seq) else 1
		
		thetype = cigar[-1]
		size = int(cigar[:-1])
		
		if thetype in ['D','H']:
			
			if size > 50:
				
				if coordinate - laststart > 500:
					
					segments.append([laststart , coordinate])
					variants.append(variant)
					
				laststart = coordinate + size
				variant = []
				
			else:
				if not ifmask:
					variant.extend([x for x in range(coordinate,coordinate + size,20)])
				
				if name in highlight_points:
					for i,seg in enumerate(highlight_points[name]):
						
						if seg[0] <= coordinate + size and seg[1] >= coordinate :
							#highlight_type[(name,seg)] = 1
							
							highlight_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
							
							
				if name in exon_points:
					for i,seg in enumerate(exon_points[name]):
						if seg[0] <= coordinate + size and seg[1] >= coordinate :
							#highlight_type[(name,seg)] = 1
							
							exon_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
							
							
			coordinate += size
			
			
		elif thetype in ['X']:
			
			
			
			if name in highlight_points:
				for i,seg in enumerate(highlight_points[name]):
					
					if seg[0] <= coordinate + size and seg[1] >= coordinate :
						#highlight_type[(name,seg)] = 1
						
						highlight_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
						
						
			if name in exon_points:
				for i,seg in enumerate(exon_points[name]):
					if seg[0] <= coordinate + size and seg[1] >= coordinate :
						#highlight_type[(name,seg)] = 1
						
						exon_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
						
			if not ifmask:
				variant.extend([x for x in range(coordinate,coordinate + size,20)])
			
			
			coordinate += size
			
		elif thetype in ['I']:
			
			if not ifmask:
				variant.append(coordinate)
				
				
			size = 1
			if name in highlight_points:
				for i,seg in enumerate(highlight_points[name]):
					
					if seg[0] <= coordinate + size and seg[1] >= coordinate :
						#highlight_type[(name,seg)] = 1
						
						highlight_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
						
						
			if name in exon_points:
				for i,seg in enumerate(exon_points[name]):
					if seg[0] <= coordinate + size and seg[1] >= coordinate :
						#highlight_type[(name,seg)] = 1
						
						exon_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
						
		else:
			coordinate += size
			
	segments.append([laststart, coordinate])
	variants.append(variant)
	
	
	return segments, variants,highlight_points,exon_points



def getused_segments(names,cigars,types):
	
	used_cigars = dict()
	for i,(name,cigar) in enumerate(zip(names,cigars)):
		
			segment_cigars = cigar.replace('<','>').split(">")[1:]
			segment_cigars = [x for x in segment_cigars if not re.match(r"[0-9H]+$", x.split(":")[1]) ]
		
			for segment_cigar in segment_cigars:
				thesize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', segment_cigar)]+[0])
				used_cigars[segment_cigar.split(":")[0].split('_')[0]] = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', segment_cigar)]+[0])
				used_cigars[segment_cigar.split(":")[0]] = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', segment_cigar)]+[0])
				
			if name in types:
					pass
			elif "_".join(name.split("_")[2:]) in types:
					pass
			else:   
					continue
		
	return used_cigars


def makemutant(inputfile, outputfile, typefile,hnames=None,hnames1=None, hnames2 = None, gff3file = "", fullfile = "" ,ifintron = 0):
	
	highlight = {}
	
	hprcmatch = "/Users/walfred/Documents/Marklab/figure4/hprcmatchresults.txt"
	matchtable = pd.read_csv(hprcmatch, sep = ',', header = None)
	matchtable = matchtable.values.tolist()
	matchresults = dict()
	for line in matchtable:
		qname = [x for x in line[1].split(";") if line[0] in x][0]
		matchresults[qname] = line[2].split(";")[0]
	matchresults ['CYP2D_group1_chr22_244'] = 'CYP2D_group1_chr22_244'
	
	usesamples = os.listdir("/Users/walfred/Documents/Marklab/figure4/results/")
	usesamples = [x.split(".")[0] for x in usesamples]
	usesamples_set = set(usesamples )
	
	
	table = pd.read_csv(inputfile, sep = '\t', header = None)
	
	names = list(table[0])	
	exons = list(table[6])
	cigars = list(table[17])
	nametocigar = {name:cigar for name, cigar in zip(names, cigars)}
	if not ifintron:
		types = {name: int(thetype.split("_")[-1]) for name, thetype, exon in zip(names, list(table[1]), exons) if exon == "Exon"}
	else:
		types = {name: int(thetype.split("_")[-1]) for name, thetype, exon in zip(names, list(table[1]), exons)}
	
	selectname = [i for i,name in enumerate(names) if name.split("_")[2] in usesamples or ("chr" in name and 'v' not in name)]
	table = table.iloc[selectname]

	names = list(table[0])	
	exons = list(table[6])
	cigars = list(table[17])
	
	used_cigars = getused_segments(names,cigars,types)

	sizes_used = []
	elements = []
	elements_ori = []
	counter = 0 
	alltypes = []

	num_hili = 0
	num_blue = 0
	num_total = 0
	for i,(name,cigar) in enumerate(zip(names,cigars)):
		
		cigar = segfiltration(cigar,used_cigars)
		print(name, cigar)
		#alltypes.append(types[name])
		
		segments, variants, highlight_points, blue_points = getsegments(cigar,highlight,exons)
		
		highlight_points = [y for seg in list(highlight_points.values()) for x in seg for y in x[2:]]
		
		blue_points =  [y for seg in list(blue_points.values()) for x in seg for y in x[2:]]
		
		variants_ori = variants
		segments_ori = segments
		#variants_ori = [[x for x in y if (x <= aldy_snprange[1] and x >= aldy_snprange[0]) or (x-9677 <= aldy_snprange[1] and x-9677 >= aldy_snprange[0]) ] for y in variants_ori ]
		
		allvariants_ori = sum(variants_ori, [])
		
		name = matchresults[name]
		
		cigar = nametocigar[name]
		
		cigar = segfiltration(cigar,used_cigars)
		
		
		alltypes.append(types[name])
		
		segments, variants, highlight_points, blue_points = getsegments(cigar,highlight,exons)
		
		#variants = [[x for x in y if (x <= aldy_snprange[1] and x >= aldy_snprange[0]) or (x-9677 <= aldy_snprange[1] and x-9677 >= aldy_snprange[0]) ] for y in variants ]
		allvariants = sum(variants, [])
		
		
		highlight_points = [[-x for x in y if x not in allvariants_ori] for y,y_ori in zip(variants, variants_ori)]
		blue_points = [[-x for x in y_ori if x not in allvariants] for y,y_ori in zip(variants, variants_ori) ]
		
		highlight_points = sum(highlight_points, [])
		blue_points =  sum(blue_points, [])
		
		elements.append([name, segments_ori, variants_ori,highlight_points, blue_points] )
		
		
		num_total += len(allvariants)
		num_blue += len(blue_points)
		num_hili += len(highlight_points)
		
		print(name,len(blue_points))
	
	
	print(num_total,num_hili,num_blue)
	
	fullsize = sum(used_cigars.values())
	#fullsize = segments[-1][-1]
	#fullsize = 1

	light_colors = [ "light"+color for color in ['yellow','skyblue' ,'gray','coral', 'seagreen' , 'steelblue', 'cyan', 'pink', 'green', 'grey','skyblue','salmon','slategray']]
	
	color_dict = {i: light_colors[i%(len(light_colors)-1) + 1] if i >0 else light_colors[0] for i in range(max(alltypes)+1)}


	# Define the height and colors of the bar segments

	plotsize = len(elements)/10

	matplotlib.rcParams["path.simplify_threshold"] = 0.00001
	matplotlib.rcParams['path.simplify'] = False

	fig,ax = plt.subplots(figsize=(2*plotsize, plotsize))

	c = 'red'
	l = 5
	L= 20


	allvariants_x = []
	allvariants_y = []

	allhighlights_x = []
	allhighlights_y = []

	allexons_x = []
	allexons_y = []
	#elements = elements[400:500]

	height = 4*plotsize/len(elements)
	width = 2*plotsize/fullsize


	theminsize = min( width,height) 
	repeattime = max( width,height) / theminsize
	repeattime = int(max(1, repeattime/10))

	highlightnames = set(hnames.split(","))
	highlightnames1 = set(hnames1.split(","))
	highlightnames2 = set(hnames2.split(","))

	labels = cl.defaultdict(list)

	typeindex = 0
	offsite = 0
	lasttype = -1
	repeattime = 1
	gapsize = 5
	print(len(elements))
	for i, (thename,segments, variants,hilis,exonpts) in enumerate(elements):
		
		if i % 100 == 0:
			print(i)
			
		thetype = alltypes[i]
		
		if thetype != lasttype:
			offsite += gapsize
			c = color_dict[typeindex]
			typeindex += 1
			
		lasttype = thetype
		
		i += offsite
		
		L = 1
		
			
		if len([x for x in  highlightnames if x in thename]):
			c = 'red'
			if len([x for x in  highlightnames1 if x in thename]):
				c = 'orange'
				
				
		elif len([x for x in  highlightnames1 if x in thename]):
			c = 'green'
		elif len([x for x in  highlightnames2 if x in thename]):
			c = 'blue'
			
			
		breaks = [x/fullsize for y in segments for x in y[:2]]
		
		
		posix = [v/fullsize for x in variants for v in x ] * repeattime
		posiy = sum([ [i+5*r*theminsize]*len(posix)  for r in range(repeattime) ], [])
		
		hilix = [-v/fullsize for v in hilis if v < 0 ] * repeattime
		hiliy = sum([ [i+5*r*theminsize]*len(hilix)  for r in range(repeattime) ], [])
		
		exonx = [-v/fullsize for v in exonpts if v < 0 ] * repeattime
		exony = sum([ [i+5*r*theminsize]*len(exonpts)  for r in range(repeattime) ], [])
		
		"""
		else:
			posix = sum([ [v/fullsize + 10*r*theminsize for x in variants for v in x]  for r in range(repeattime) ], [])
			posiy =  [i+r*theminsize]* len(posix) * repeattime
			
			hilix = sum([ [-v/fullsize + 10*r*theminsize for v in hilis if v < 0 ]  for r in range(repeattime) ], [])
			hiliy = [i+r*theminsize]* len(posix) * repeattime
			
			exonx = sum([ [-v/fullsize + 10*r*theminsize for v in exonpts if v < 0 ]  for r in range(repeattime) ], [])
			exony = [i+r*theminsize]* len(posix) * repeattime
		"""	
		
		allvariants_x.extend(posix)
		allvariants_y.extend( posiy)
		
		allhighlights_x.extend(hilix)
		allhighlights_y.extend(hiliy)
		
		
		allexons_x.extend(exonx)
		allexons_y.extend(exony)
		
		last_coordinate = 0
		for index, coordinate in enumerate(breaks):
			
			if index % 2 == 0:
				l = 0.3
				ax.plot([last_coordinate, coordinate], [i, i], color=c, linewidth=l, linestyle='-.')
			else:
				l = 3
				ax.plot([last_coordinate, coordinate], [i, i], color=c, linewidth=l, linestyle='-')
				
			last_coordinate = coordinate
			
		labels[thetype].append(i)
		
		
	custom_grey = (0.100, 0.100, 0.100)  # RGB values between 0 and 1
	ax.scatter(allvariants_x, allvariants_y, marker='.',color =  'grey' ,edgecolor='grey',s = 3,zorder=10,linestyle='None',alpha=1, cmap='viridis')
	ax.scatter(allexons_x, allexons_y, marker='.',color =  'blue'  ,edgecolor='blue',s = 30,zorder=10,linestyle='None',alpha=1, cmap='viridis')
	ax.scatter(allhighlights_x, allhighlights_y, marker='.',color =  'red'  ,edgecolor='red',s = 30,zorder=10,linestyle='None',alpha=1, cmap='viridis')
	# Show the plot

	plt.axis('off')
	
	if len(gff3file):
		
		if len(fullfile) == 0:
			fullfile = inputfile
			
		plotexons(fig, ax, i , fullfile, gff3file, used_cigars)
		
	plt.savefig(outputfile)	
	
	
	
	
	
def main(args):
	
		makemutant(args.input, args.output, args.type, args.name,args.name1, args.name2, args.gff, args.full, args.intron)
	
def run():
		"""
				Parse arguments and run
		"""
		parser = argparse.ArgumentParser(description="program determine psuedogene")
		parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
		parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
		parser.add_argument("-t", "--type", help="path to output file", dest="type",type=str, default = "")
		parser.add_argument("-n", "--name", help="path to output file", dest="name",type=str, default = "highlight")
		parser.add_argument("-n2", "--name2", help="path to output file", dest="name2",type=str, default = "highlight")
		parser.add_argument("-n1", "--name1", help="path to output file", dest="name1",type=str, default = "highlight")
		parser.add_argument("-g", "--gff", help="path to output file", dest="gff",type=str, default = "")
		parser.add_argument("-f", "--full", help="path to output file", dest="full",type=str, default = "")
		parser.add_argument("-intron", "--intron", help="path to output file", dest="intron",type=int, default = 0)
	
		parser.set_defaults(func=main)
		args = parser.parse_args()
		args.func(args)
	
	
if __name__ == "__main__":
		run()
	
	
