#!/usr/bin/env python3

import re
import numpy as np
import collections as cl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

class UPGMANode:
	
	def __init__(self, left=None, right=None, up_dist=0.0, down_dist=0.0):
		
		self.left = left
		self.right = right
		self.up_dist = up_dist
		self.down_dist = down_dist
		self.nleaves = 1
		
		if type(left) is type(UPGMANode):
			self.nleaves += left.nleaves + right.nleaves
		else:
			self.nleaves = 1
			
			
	def leaves(self) -> list:
		
		if self is None:
			return []
		elif self.right == None:
			return [self.left]
		
		return self.left.leaves() + self.right.leaves()	
	
	def to_newick(self) -> str:
		
		if self.right == None:
			return self.left + ":" + str(self.up_dist)
		else:
			return ( "(" + ",".join([x.to_newick() for x in [self.left, self.right]]) + "):"+ str(self.up_dist) )
		
class UPGMA:
	
	
	def __init__(self, dist_matrix: np.ndarray, header: list):
		
		self.distances = dist_matrix
		self.header = header
		self.exclude = [0 for x in range(len(header))]+[1 for x in range(len(header))]
		self.build_tree(self.distances, self.header)
		
	def getmindist(self)-> tuple:
		
		min_row, min_col = 0,0
		
		curr_min = np.inf
		for row,values in enumerate(self.work_matrix):
			
			if self.exclude[row]:
				continue
			
			for col, value in enumerate(values):
				
				if self.exclude[col]:
					continue
				
				if value == 0:
					
					return (row, col)
				
				elif value < curr_min:
					
					curr_min = value
					min_row = row
					min_col = col
					
					
		return (min_row, min_col)
	
	
	def build_tree(self, dist_matrix: np.ndarray, header: list) -> UPGMANode:
		
		nodes = [UPGMANode(taxon) for taxon in header]
		
		headertoindex = {name:i for i,name in enumerate(header)}
		
		self.size = 2*len(header)
		
		self.work_matrix = np.array([row+[np.inf]*len(row) for row in dist_matrix.tolist()]+[[np.inf]*(2*len(row)) for row in dist_matrix.tolist()], dtype=float)
		np.fill_diagonal(self.work_matrix, np.inf)
		
		
		for turn in range(len(nodes)-1):
			
			least_id = self.getmindist()
			least_dist = self.work_matrix[least_id[0], least_id[1]]
			node1, node2 = nodes[least_id[0]], nodes[least_id[1]]
			self.exclude[least_id[0]] = 1
			self.exclude[least_id[1]] = 1
			
			new_node = UPGMANode(node2, node1)
			nodes.append(new_node)
			self.exclude[len(nodes)-1] = 0
			node1.up_dist = least_dist / 2 - node1.down_dist
			node2.up_dist = least_dist / 2 - node2.down_dist
			new_node.down_dist = least_dist / 2
			
			# create new working distance matrix
			self.update_distance( nodes, least_id)
			
		self.tree = new_node 
		
	def update_distance(self, nodes: list, least_id: tuple) -> np.ndarray:
		
		length = len(nodes)
		nleaves1, nleaves2 = nodes[least_id[0]].nleaves, nodes[least_id[1]].nleaves
		nleaves  = nleaves1 + nleaves2 
		
		for i in range(length-1):
			
			if self.exclude[i]:
				continue
			
			
			self.work_matrix[i,length-1] = ( self.work_matrix[i][least_id[0]]*nleaves1 + self.work_matrix[i][least_id[1]]*nleaves2 ) / nleaves
			self.work_matrix[length-1,i] = self.work_matrix[i,length-1]
			
			

def ordermatch(l1,l2):
	
	l1 = [0] + l1 + [0]
	l2 = [0] + l2 + [0]
	
	l1_path = [(l1[i-1],l1[i]) for i in range(1,len(l1))]
	l2_path = [(l2[i-1],l2[i]) for i in range(1,len(l2))]
	
	match = 0 
	for path in l1_path:
		if path in l2_path:
			match +=1 
			l2_path.remove(path)
	
	return 1.0 - (2*match)/(len(l1_path) + len(l2_path))
	
	

def projectiontree(dm,l):
	
	header = ["_"+str(x)+"_" for x in range(0,l)]
	
	tree = UPGMA(dm, header).tree
	
	
	return tree.to_newick()

svtypes = [[5,6,7,8,9,10],[11], [12,13,14,15], [16],[17,18,19,20],[21,22,23,24], [25,26,27,28,29,30,31],[32],[33],[34]]

thefile = "./AMY_group1_AMY1AOOOAMY1BOOOAMY1C.fasta_annotate.fa_annotatesummary.txt"

with open(thefile, mode = 'r') as f:
	text = [line.strip().split('\t') for line in f.readlines()]
	

singleton_sv_types =set()
typetosvtype = dict()
for i,thetype in enumerate(svtypes):
	for index in thetype:
		typetosvtype[index] = i+1
	if len(thetype) == 1:
		singleton_sv_types.add(thetype[0])


haplos = cl.defaultdict(list)
for line in text:
	
	cladeindex = int(line[1].split("_")[-1])
	cladeindex = typetosvtype.get(cladeindex, 0)
	
	original = line[7]
	haplo, region, strd = original.split(":")[0],original.split(":")[1][:-1], original.split(":")[1][-1]
	
	start,end = region.split("-")
	
	haplos[haplo].append([start, end,cladeindex, strd])

haplo_types = []
for haplo, data in haplos.items():
	
	strds = [x[-1] for x in data]
	numposi = sum(1 for x in strds if x=="+")
	if 2*numposi > len(strds):
		thestrd = 1
	else:
		thestrd = -1
	
	data = [x for x in data if abs(x[2]) > 0]
	data = sorted(data, key = lambda x:int(x[0])*strd)
	if thestrd == -1:
		for i in range(len(data)):
			data[i][-1] = '+' if data[i][-1] == '-' else '-'
			
	data = [x[2] if x[-1] == '+' else -x[2] for x in data]
	if haplo in ["chr1", "NC_060925.1"]:
		print(haplo, data)
		
	if len(data):
		haplo_types.append(data)
	
matrix = np.ones((len(haplo_types),len(haplo_types))) 

for i, rowi in enumerate(haplo_types):
	
	matrix[i, i] = 1.0
	for j, rowj in enumerate(haplo_types[(i+1):]):
		
		matrix[i, i+j+1] = ordermatch(rowi,rowj)
		matrix[i+j+1, i] = matrix[i, i+j+1]
	

thetree = projectiontree(matrix,len(haplo_types))

allorders = re.findall(r'_\d+_',thetree)
allorders = [int(x[1:-1]) for x in allorders]


haplo_types_order = [haplo_types[i] for i in allorders ]

type_count = cl.defaultdict(int)
for haplo in haplo_types_order:
	
	type_count[tuple(haplo)]+=1

type_sort = sorted(list(type_count.keys()), key = lambda x: len(x))

haplo_types_order = []
for haplo in type_sort :
	for x in range(type_count[haplo]):
		haplo_types_order.append(list(haplo))
	
data = haplo_types_order


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

fig, ax = plt.subplots(figsize=(20, 50))
lastrow = None
i = 0  # Initialize row counter for spacing between rows

for index, row in enumerate(data):
	
	lastrow = row
	i += 1
	
	for j, value in enumerate(row):
		
		color = distinct_colors[abs(value)-1]  # Determine color based on value
		start_x = j
		end_x = start_x + 1 # Adjust length and direction
		
		if value < 0:
			temp = start_x
			start_x = end_x
			end_x = temp
			
			
		arrow = patches.FancyArrowPatch(
			(2*start_x, i*5), (2*end_x, i*5),
			mutation_scale=20,
			#arrowstyle='->,head_width=0.3,head_length=1',  # Change arrow style and head size
			facecolor=color,  # Arrow fill color
			linewidth=0.00,  # Frame thickness
			#edgecolor='black',  # Black frame around the arrow
		)
		ax.add_patch(arrow)  # Add arrow to the axes
		
		
# Remove axis and labels
ax.axis('off')  # Turn off the axis

# Set limits to fit the arrows properly
ax.set_xlim(-2, 2 * len(max(data, key=len))+4)
ax.set_ylim(-5, 5*i+5)

plt.savefig('./AMY_chrom.pdf')
