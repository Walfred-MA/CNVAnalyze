#!/usr/bin/env python3

# -i AMY_group1_AMY1AOOOAMY1BOOOAMY1C.fasta_annotate.fa_annotatesummary_AMY1.txt -g AMY.gff3 -f AMY_group1_AMY1AOOOAMY1BOOOAMY1C.fasta_annotate.fa_annotatesummary.txt -o AMY_major.png

import re
import collections as cl
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
			
				if thesize >0 and thesize < 0:
					
						fullsize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', thecigar)]+[0])
					
						thecigar = "{}H".format(fullsize )
					
				newcigars.append(thename+":"+thecigar)
			
		newcigars= "".join(newcigars)
	
		return newcigars

def getgenelocation(table, gff3file, used_cigars, segment_orders):
	
	segment_exonalign = cl.defaultdict(lambda : cl.defaultdict(list))
	gffrecord = cl.defaultdict(lambda : cl.defaultdict(list))
	if len(gff3file):
		
		with open(gff3file, mode = 'r') as f:
			for line in f:
				line = line.strip().split()
				
				if line[2] == 'exon':
					genename = [x for x in line[-1].split(";") if x.startswith("gene_name=")][0][10:]
					if genename != "AMY2A":
						gffrecord[line[0]][genename].append([int(line[3]), int(line[4])])
		
		gffrecords = {thechr:{gene:overlapexons(coordi) for gene, coordi in record.items()} for thechr, record in gffrecord.items()}
		
		for row in table.values.tolist():
			
			locus = row[7]
			refchr = locus.split(":")[0]
			if refchr not in gffrecords:
				continue
			
			strd = locus[-1]
			
			thepath = row[-3]
			current_order = [x[1:] for x in re.findall('[<>][0-9_]+',thepath)]
			
			cigar = row[-2]
			segment_cigars = re.findall('[<>][0-9_A-Za-z:]+',segfiltration(cigar))
			
			segment_cigars = {x.split(":")[0][1:]:x for x in segment_cigars}
			
			segment_cigars = [segment_cigars[x] for x in current_order]
			
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
				
				
				
				if strd == '+':
					exons = sum([ [ max(0, x[0] - start), x[1] - start ] for x in exons],[])
				else:
					exons = sum([ [ max(0, end - 1 - x[1]), end - x[0] ] for x in exons[::-1]],[])
				
				exons_onsegments = maptoseg(exons, segment_qranges)
				
				print(gene,exons_onsegments)
				
				for cigar, exons_onsegment, size in zip(segment_cigars, exons_onsegments, segment_alignsizes):
					
					if len(exons_onsegment) == 0:
						continue
					
					name, cigar = cigar.split(":")
					
					
					if name[1:] not in used_cigars:
						continue
	
					if name[0] == '<':
						exons_onsegment = sorted([size[1]-x for x in exons_onsegment])
					
					
					exons_onsegment_rposi = cigar_findrposi(cigar, exons_onsegment)
				
									
					segment_exonalign[gene][name[1:]].extend(exons_onsegment_rposi)
			
		for gene, segment_exonalign_ in segment_exonalign.items():	
		
			for name, exonalign in segment_exonalign_.items():
				
				exonalign = overlapexons(exonalign)
				
				
				exonalign = [x for x in exonalign if x[1]-x[0] > 30]
			
			
				segment_exonalign[gene][name] = exonalign
	
	return segment_exonalign
	

def plotexons(fig, ax, lelement, inputfile, gff3file, used_cigars, segment_orders):
	
	table = pd.read_csv(inputfile, sep = '\t', header = None)
	
	exonlocations = getgenelocation(table, gff3file, used_cigars, segment_orders)
	
	exonlocations_order = {gene: [exonlocation[segment] for segment in segment_orders ]  for gene,exonlocation in exonlocations.items() }
	
	print(exonlocations_order)
	c = 'red'
	l = 2
	L= 20
	chromy = lelement + 10
	fullsize = sum(list(used_cigars.values()))
	for gene, exonlocation in exonlocations_order.items():
		
		if len(sum(exonlocation, [])) == 0:
			continue

		
		name_add = cl.defaultdict(int)
		leftend = 0
		exons_posi = []
		
		repeatindex = set()
		for name, pathexon in zip(segment_orders, exonlocation):
			
			if len(pathexon) :
				name_add[name] += 1	
				
			size = used_cigars[name]
			rightend = leftend + size
			
			for exon in pathexon:
				exons_posi.append([ (leftend + exon[0]) / fullsize, (leftend + exon[1]) / fullsize])
			
			if name_add[name] > 1:
				
				repeatindex.add(len(exons_posi) -1)
				
			leftend = rightend
		
		exons_posi_ = sum([x for i,x in enumerate(exons_posi) if i not in repeatindex], [])
		exons_posi_min, exons_posi_max = min(exons_posi_), max(exons_posi_)
		ax.plot([exons_posi_min, exons_posi_max], [chromy, chromy], color=c, linewidth=l, linestyle='-')
		for i,exon in enumerate(exons_posi):
			if i not in repeatindex:
				ax.plot(exon, [chromy, chromy], color=c, linewidth=L, linestyle='-')
		
		if len(repeatindex):
			for i,exon in enumerate(exons_posi):
				if i in repeatindex:
					ax.plot(exon, [lelement + 10+7*len(exonlocations_order), lelement + 10+7*len(exonlocations_order)], color=c, linewidth=L, linestyle='-')
				
	
		chromy += 7
	
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

def getsegments(cigarstr,highlight):
	
		cigars = re.findall(r'[><0-9_:]*\d+[HMX=ID]', cigarstr)
	
		variants = []
		segments = []
		highlight_points = {}
		highlight_coordinate = -1
		highlight_type = cl.defaultdict(int)
	
		oldname = ""
		variant = []
		coordinate = 0
		laststart = 0
		findposition = -1
		for cigar in  cigars:
			
				if ":" in cigar:
					
						name,cigar = cigar.split(":")
					
						name = name[1:]
						if name in highlight:
							
								highlight_points[name] = highlight[name] + coordinate
								highlight_coordinate = highlight[name] + coordinate
							
							
				thetype = cigar[-1]
				size = int(cigar[:-1])
			
				if thetype in ['D','H']:
					
						if size > 0:
							
								if coordinate - laststart > 300:
									
										segments.append([laststart , coordinate])
										variants.append(variant)
									
								laststart = coordinate + size
								variant = []
						else:
							variant.extend([x for x in range(coordinate,coordinate + size,20)])
						
						if highlight_coordinate >= coordinate and highlight_coordinate < coordinate + size:
								highlight_type[(name,highlight_coordinate)] = -1
							
						coordinate += size
					
					
				elif thetype in ['X']:
					
						if highlight_coordinate >= coordinate and highlight_coordinate < coordinate + size:
							
								highlight_type[(name,highlight_coordinate)] = 1
								highlight_points[name] *= -1
							
							
						variant.extend([x for x in range(coordinate,coordinate + size,20)])
						coordinate += size
					
				elif thetype in ['I']:
					
						variant.append(coordinate)
				else:
						coordinate += size
		
		
		segments.append([laststart, coordinate])
		variants.append(variant)
	
	
		return segments, variants,highlight_points



def makemutant(inputfile, outputfile, typefile,hnames=None,hnames1=None, hnames2 = None, gff3file = "", fullfile = "" ,ifintron = 0):
	
		svtypes = [[5,6,7,8,9,10],[11], [12,13,14,15], [16],[17,18,19,20],[21,22,23,24], [25,26,27,28,29,30,31],[32],[33],[34]]
	
		singleton_sv_types =set()
		typetosvtype = dict()
		for i,thetype in enumerate(svtypes):
			for index in thetype:
				typetosvtype[index] = i+1
			if len(thetype) == 1:
				singleton_sv_types.add(thetype[0])
	
		highlight = {} 
		
		table = pd.read_csv(inputfile, sep = '\t', header = None)
		
		names = list(table[0])		
		exons = list(table[6])
		
		if not ifintron:
			types = {name: int(thetype.split("_")[-1]) for name, thetype, exon in zip(names, list(table[1]), exons) if exon == "Exon"}
		else:
			types = {name: int(thetype.split("_")[-1]) for name, thetype, exon in zip(names, list(table[1]), exons)}
		
		types = {name: typetosvtype[thetype] for name, thetype in types.items()}
		
		cigars = list(table[17])
	
	
		used_cigars = dict()
		for i,(name,cigar) in enumerate(zip(names,cigars)):
			
				segment_cigars = cigar.replace('<','>').split(">")[1:]
				segment_cigars = [x for x in segment_cigars if not re.match(r"[0-9H]+$", x.split(":")[1]) ]
				
				for segment_cigar in segment_cigars:
					
					thesize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', segment_cigar)]+[0])
					
					if thesize>= 300:
						used_cigars[segment_cigar.split(":")[0].split('_')[0]] = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', segment_cigar)]+[0])
						used_cigars[segment_cigar.split(":")[0]] = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', segment_cigar)]+[0])
			
				if name in types:
						pass
				elif "_".join(name.split("_")[2:]) in types:
						pass
				else:   
						continue

		cigars0 = cigars[0]
		segment_cigars = re.findall('[<>][0-9_A-Za-z:]+',cigars0)
	
		segment_orders = []
		for segment_cigar in segment_cigars:
			segment_orders.append(segment_cigar.split(':')[0][1:])
			
		segment_orders = [x.split('_')[0] for x in segment_orders if x in used_cigars]

		
		sizes_used = []
		elements = []
		counter = 0 
		alltypes = []
		for i,(name,cigar) in enumerate(zip(names,cigars)):
			
				cigar = segfiltration(cigar,used_cigars)
				if name in types:
						alltypes.append(types[name])
				elif "_".join(name.split("_")[2:]) in types:
						alltypes.append(types["_".join(name.split("_")[2:])])
				else:
						print("skipping: ",name)
						continue
						#alltypes.append(0)
			
				segments, variants, highlight_points = getsegments(cigar,highlight)
			
				if len(highlight_points):
						counter += 1
					
				elements.append([names[i], segments, variants,list(highlight_points.values())])
			
		#fullsize = segments[-1][-1]
		fullsize = sum(used_cigars.values())
		
		#light_colors = [ "light"+color for color in ['yellow','skyblue' ,'gray','coral', 'seagreen' , 'steelblue', 'cyan', 'pink', 'green', 'grey','skyblue','salmon','slategray']]
	
		
	
		distinct_colors = [
			'#e6194b',  # Red
			'#3cb44b',  # Green
			'#ffe119',  # Yellow
			'#4363d8',  # Blue
			"pink",
			'#f58231',  # Orange
			'#911eb4',  # Purple
			'#46f0f0',  # Cyan
			'#f032e6',  # Magenta
			'#bcf60c'   # Lime
		]
		light_colors = [ color for color in distinct_colors]
		
		color_dict = {i: light_colors[(i-1)%len(light_colors)] for i in range(max(alltypes)+1)}
	
	
		# Define the height and colors of the bar segments
		offsite = 0
		gapsize = 10
		lasttype = -1
		
		
		plotsize = min(200,len(elements)/10)
		fig,ax = plt.subplots(figsize=(plotsize, plotsize))
		
		test = 0
		if len(gff3file) and test == 1:
			
			if len(fullfile) == 0:
				fullfile = inputfile
				
			plotexons(fig, ax, i , fullfile, gff3file, used_cigars)
		
	
		allvariants_x = []
		allvariants_y = []
	
		allhighlights_x = []
		allhighlights_y = []
		#elements = elements[400:500]
	
		height = 4*plotsize/len(elements)
		width = 2*plotsize/fullsize
	
	
		theminsize = min( width,height) 
		repeattime = max( width,height) / theminsize
		repeattime = int(max(1, repeattime/10))
	
		highlightnames = set(hnames.split(","))
		highlightnames1 = set(hnames1.split(","))
		highlightnames2 = set(hnames2.split(","))
	
		
		lasttype = -1
		repeattime = 1
		for i, (thename,segments, variants,hilis) in enumerate(elements):
			
				if i % 100 == 0:
						print(i)
					
					
				thetype = alltypes[i]
				c = color_dict[thetype]
				if c != lasttype:
						offsite += gapsize
					
				lasttype = c
			
				i += offsite
			
				L = 1
				if "chr" in thename and "v" not in thename :
						#c = c.replace("light","")
						c = adjust_lightness(c, amount=0.7)
						c = "darkblue"
						L = 3
				#elif "NC_0609" in thename:
						#c = 'red'
					
					
				if len([x for x in  highlightnames if x in thename]):
						c = 'red'
						if len([x for x in  highlightnames1 if x in thename]):
								c = 'orange'
							
				elif len([x for x in  highlightnames1 if x in thename]):
						c = 'green'
				elif len([x for x in  highlightnames2 if x in thename]):
						c = 'blue'
					
					
				breaks = [x/fullsize for y in segments for x in y[:2]]
			
			
				posix = [v/fullsize for x in variants for v in x] * repeattime
				posiy = sum([ [i+10*r*theminsize]*len(posix)  for r in range(repeattime) ], [])
			
				hilix = [-v/fullsize for v in hilis if v < 0 ] * repeattime
				hiliy = sum([ [i+10*r*theminsize]*len(hilix)  for r in range(repeattime) ], [])
			
				"""
				else:
						posix = sum([ [v/fullsize + 10*r*theminsize for x in variants for v in x]  for r in range(repeattime) ], [])
						posiy =  [i+r*theminsize]* len(posix) * repeattime

						hilix = sum([ [-v/fullsize + 10*r*theminsize for v in hilis if v < 0 ]  for r in range(repeattime) ], [])
						hiliy = [i+r*theminsize]* len(posix) * repeattime
				"""
			
			
				allvariants_x.extend(posix)
				allvariants_y.extend( posiy)
			
				allhighlights_x.extend(hilix)
				allhighlights_y.extend(hiliy)
			
				last_coordinate = 0
				for index, coordinate in enumerate(breaks):
					
						if index % 2 == 0:
								l = 0.2
								ax.plot([last_coordinate, coordinate], [i, i], color=c, linewidth=l, linestyle='-.')
						else:
								l = 2
								ax.plot([last_coordinate, coordinate], [i, i], color=c, linewidth=l, linestyle='-')
							
						last_coordinate = coordinate
					
					
		ax.scatter(allvariants_x, allvariants_y, marker='.',color =  'black'  ,s = 5,zorder=3,linestyle='None')
		ax.scatter(allhighlights_x, allhighlights_y, marker='.',color =  'red'  ,s = 10,zorder=3,linestyle='None')
		# Show the plot
		if len(gff3file):
			
			if len(fullfile) == 0:
				fullfile = inputfile
				
			plotexons(fig, ax, i , fullfile, gff3file, used_cigars, segment_orders)
	
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
	
