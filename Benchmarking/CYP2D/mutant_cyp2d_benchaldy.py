#!/usr/bin/env python3

#/Users/walfred/anaconda3/bin/python3 $filename -i /Users/walfred/Documents/Marklab/figure4/CYP2D_group1_CYP2D6OOOCYP2D7OOOCYP2D8P.fasta_annotate.fa_lineargraph.gaf -t /Users/walfred/Documents/Marklab/figure4/CYP2D_group1_CYP2D6OOOCYP2D7OOOCYP2D8P.fasta_annotate.fa_types.txt -o ./CYP2D6_bench.png

import re
import collections as cl
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
import os 
import argparse

def sizefiltration(cigar):
	
	cigars = re.findall('[<>][0-9_A-Za-z:]+',cigar)
	
	newcigars = []
	for cigarstr in cigars:
		
		thename,thecigar = cigarstr.split(":")
		
		thesize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=]', thecigar)]+[0])
		
		if thesize >0 and thesize < 200:
			
			fullsize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', thecigar)]+[0])
			
			thecigar = "{}H".format(fullsize )
			
		newcigars.append(thename+":"+thecigar)
		
	newcigars= "".join(newcigars)
	
	return newcigars

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

def getsegments(cigarstr,highlight,exons,pathseqs):
	
	cigars = re.findall(r'[><0-9_:]*\d+[HMX=ID]', cigarstr)
	
	
	variants = []
	segments = []
	highlight_points = {}
	blue_points = {}
	
	highlight_coordinate = -1
	highlight_type = cl.defaultdict(int)
	
	oldname = ""
	variant = []
	
	localpath_offsite = 0
	coordinate = 0
	laststart = 0
	findposition = -1
	currentseq = ""
	
	for cigar in  cigars:
		
		if ":" in cigar:
			
			localpath_offsite = coordinate
			
			name,cigar = cigar.split(":")
			
			currentseq = pathseqs[name[1:].split("_")[0]]
			
			name = name[1:]
			
			if name in highlight:
				
				highlight_points[name] = [[x + coordinate for x in seg] for seg in highlight[name] ]
			
			if name in exons:
				
				blue_points[name] =  [[x + coordinate for x in seg] for seg in exons[name] ]
				
		if len(currentseq) == 0:
			continue
		
		
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
				variant.extend([x for x in range(coordinate,coordinate + size,20)])
				if name in highlight_points:
					for i,seg in enumerate(highlight_points[name]):
						
						if seg[0] <= coordinate + size and seg[1] >= coordinate :
							#highlight_type[(name,seg)] = 1
							
							highlight_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
							
							
				if name in blue_points:
					for i,seg in enumerate(blue_points[name]):
						if seg[0] <= coordinate + size and seg[1] >= coordinate :
							#highlight_type[(name,seg)] = 1
							
							blue_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
			
			"""
			for seg in highlight_points:
				if highlight_coordinate[0] >= coordinate and highlight_coordinate[1] <= coordinate + size:
					highlight_type[(name,highlight_coordinate)] = -1
					break
			"""
					
				
			coordinate += size
			
			
		elif thetype in ['X']:
			
			
			
			if name in highlight_points:
				for i,seg in enumerate(highlight_points[name]):
					
					if seg[0] <= coordinate + size and seg[1] >= coordinate :
						#highlight_type[(name,seg)] = 1
					
						highlight_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
					
			
			if name in blue_points:
				for i,seg in enumerate(blue_points[name]):
					if seg[0] <= coordinate + size and seg[1] >= coordinate :
						#highlight_type[(name,seg)] = 1
						
						blue_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
						

			variant.extend([x for x in range(coordinate,coordinate + size,20) if currentseq[coordinate - localpath_offsite].isupper()])
			
			
			coordinate += size
			
		elif thetype in ['I']:
			
			if currentseq[coordinate - localpath_offsite].isupper():
				variant.append(coordinate)
			
			
			size = 1
			if name in highlight_points:
				for i,seg in enumerate(highlight_points[name]):
					
					if seg[0] <= coordinate + size and seg[1] >= coordinate :
						#highlight_type[(name,seg)] = 1
						
						highlight_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
					
				
			if name in blue_points:
				for i,seg in enumerate(blue_points[name]):
					if seg[0] <= coordinate + size and seg[1] >= coordinate :
						#highlight_type[(name,seg)] = 1
						
						blue_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
						
			
		else:
			coordinate += size
			
	segments.append([laststart, coordinate])
	variants.append(variant)
	
	
	return segments, variants,highlight_points,blue_points



def makemutant(inputfile, outputfile, typefile,hnames=None,hnames1=None, hnames2 = None):
	
		
	highlight= {}
	exons = {}
	pathsizes = {}
	pathseqs = dict()
	lines = []
	with open(inputfile, mode = 'r') as f:
		
		for line in f:
		
			if len(line.strip()) == 0:
				continue
			
			if line[0] == 'L':
				lines.append(line.strip().split('\t'))
			elif line[0] == 'S'	:
				line = line.strip().split()
				name = line[1][2:]
				pathsizes[name] = int(line[-2].split(':')[-1])
				
				if "chr" in line[-1] and name in ['7','8','9']:
					pathseqs[name] = line[2]
				else:
					pathseqs[name] = ""
		
	
	table = pd.DataFrame.from_records(lines)
	
	types = cl.defaultdict(int)
	if len(typefile):
		with open(typefile, mode = 'r') as f:
			text = [x.split(",") for x in f.read().splitlines() if len(x)]
			for i, names in enumerate(text):
				for name in names:
					types[name] = i + 1
	
	hprcmatch = "/Users/walfred/Documents/Marklab/figure4/hprcmatchresults.txt"
	matchtable = pd.read_csv(hprcmatch, sep = ',', header = None)
	
	matchtable = matchtable.values.tolist()

	matchresults = dict()
	for line in matchtable:
		
		qname = [x for x in line[1].split(";") if line[0] in x][0]
		matchresults[qname] = line[2].split(";")[0]
	
	matchresults ['CYP2D_group1_chr22_244'] = 'CYP2D_group1_chr22_244'
	
	aldy_snprange= [3271, 9336]
	names = list(table[1])
	cigars = list(table[5])
	nametocigar = {name:cigar for name, cigar in zip(names, cigars)}
	sizes = [int(x.split(':')[-1]) for x in list(table[6])]
	
	sizes_used = []
	elements = []
	elements_ori = []
	counter = 0 
	alltypes = []
	
	usesamples = os.listdir("/Users/walfred/Documents/Marklab/figure4/results/")
	
	num_hili = 0
	num_blue = 0
	num_total = 0
	for i,(name,cigar) in enumerate(zip(names,cigars)):
		
		cigar = sizefiltration(cigar)
		
		
		if name.split("_")[2]+".final.cram.out" not in usesamples and ("chr" not in name or 'v' in name):
			continue
		
		if ">7" not in cigar and ">8" not in cigar and  ">9" not in cigar :
			continue
		
		
		if name in types:
			alltypes.append(types[name])
		elif "_".join(name.split("_")[2:]) in types:
			alltypes.append(types["_".join(name.split("_")[2:])])
		else:
			pass
			#continue
			#alltypes.append(0)
		print(name)
		
		segments, variants, highlight_points, blue_points = getsegments(cigar,highlight,exons,pathseqs)
		
		highlight_points = [y for seg in list(highlight_points.values()) for x in seg for y in x[2:]]
				
		blue_points =  [y for seg in list(blue_points.values()) for x in seg for y in x[2:]]
			
		variants_ori = variants
		
		variants_ori = [[x for x in y if (x <= aldy_snprange[1] and x >= aldy_snprange[0]) or (x-9677 <= aldy_snprange[1] and x-9677 >= aldy_snprange[0]) ] for y in variants_ori ]
		
		allvariants_ori = sum(variants_ori, [])
		
		
		
		name = matchresults[name]
		
		cigar = nametocigar[name]
		
		cigar = sizefiltration(cigar)
		
		
		if name in types:
			alltypes.append(types[name])
		elif "_".join(name.split("_")[2:]) in types:
			alltypes.append(types["_".join(name.split("_")[2:])])
		else:
			pass
			#continue
			#alltypes.append(0)
			
		segments, variants, highlight_points, blue_points = getsegments(cigar,highlight,exons,pathseqs)
		
		variants = [[x for x in y if (x <= aldy_snprange[1] and x >= aldy_snprange[0]) or (x-9677 <= aldy_snprange[1] and x-9677 >= aldy_snprange[0]) ] for y in variants ]
		allvariants = sum(variants, [])
		
		
		highlight_points = [[-x for x in y if x not in allvariants_ori] for y,y_ori in zip(variants, variants_ori)]
		blue_points = [[-x for x in y_ori if x not in allvariants] for y,y_ori in zip(variants, variants_ori) ]
		
		highlight_points = sum(highlight_points, [])
		blue_points =  sum(blue_points, [])
		
		elements.append([name, segments, variants,highlight_points, blue_points] )
		
		
		num_total += len(allvariants)
		num_blue += len(blue_points)
		num_hili += len(highlight_points)
		
		print(name,len(blue_points))
		
	print(num_total,num_hili,num_blue)
	fullsize = segments[-1][-1]
	#fullsize = 1
	
	light_colors = [ "light"+color for color in ['white','yellow','skyblue' ,'gray','coral', 'seagreen' , 'steelblue', 'cyan', 'pink', 'green', 'grey','skyblue','salmon','slategray']]
	
	color_dict = {i: light_colors[i%(len(light_colors)-1) + 1] if i >0 else light_colors[0] for i in range(max(alltypes)+1)}
	
	
	# Define the height and colors of the bar segments

	plotsize = len(elements)/10
	
	matplotlib.rcParams["path.simplify_threshold"] = 0.00001
	matplotlib.rcParams['path.simplify'] = False
	
	fig,ax = plt.subplots(figsize=(plotsize, plotsize))
	
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
	
	offsite = 0
	lasttype = -1
	repeattime = 1
	gapsize = 0
	for i, (thename,segments, variants,hilis,exonpts) in enumerate(elements):
		
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
			c = c.replace("light","")
			L = 3
			
		
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
		posiy = sum([ [i+10*r*theminsize]*len(posix)  for r in range(repeattime) ], [])
		
		hilix = [-v/fullsize for v in hilis if v < 0 ] * repeattime
		hiliy = sum([ [i+10*r*theminsize]*len(hilix)  for r in range(repeattime) ], [])
		
		exonx = [-v/fullsize for v in exonpts if v < 0 ] * repeattime
		exony = sum([ [i+10*r*theminsize]*len(exonpts)  for r in range(repeattime) ], [])
			
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
				l = L*0.2
				ax.plot([last_coordinate, coordinate], [i, i], color=c, linewidth=l, linestyle='-.')
			else:
				l = L*2
				ax.plot([last_coordinate, coordinate], [i, i], color=c, linewidth=l, linestyle='-')
				
			last_coordinate = coordinate
		
		labels[thetype].append(i)
	

		
	ax.scatter(allvariants_x, allvariants_y, marker='.',color =  'grey'  ,edgecolor='grey',s = 10,zorder=10,linestyle='None',alpha=1, cmap='viridis')
	ax.scatter(allexons_x, allexons_y, marker='.',color =  'blue'  ,edgecolor='blue',s = 10,zorder=10,linestyle='None',alpha=1, cmap='viridis')
	ax.scatter(allhighlights_x, allhighlights_y, marker='.',color =  'red'  ,edgecolor='red',s = 10,zorder=10,linestyle='None',alpha=1, cmap='viridis')
	# Show the plot
	
	
	plt.axis('off')
	plt.savefig(outputfile)
	
	
def main(args):
	
	makemutant(args.input, args.output, args.types, args.name,args.name1, args.name2)
	
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine psuedogene")
	parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
	parser.add_argument("-t", "--types", help="path to output file", dest="types",type=str, default = "")
	parser.add_argument("-n", "--name", help="path to output file", dest="name",type=str, default = "highlight")
	parser.add_argument("-n2", "--name2", help="path to output file", dest="name2",type=str, default = "highlight")
	parser.add_argument("-n1", "--name1", help="path to output file", dest="name1",type=str, default = "highlight")
	
	
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
	
	
