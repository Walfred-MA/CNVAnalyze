#!/usr/bin/env python3


import pandas as pd
import os 
import argparse
import sys
import numpy as np
import multiprocessing as mul
import math
import collections as cl


def withinvar(matrix, indexes):

        avars = []
        rvars = []
        for index, i in enumerate(indexes):

                size1 = max(1.0 , matrix[i,i])
                #row = dm[ i]

                for j in indexes[:index]:

                        size2 = max(1.0, matrix[j,j])
                        match = matrix[i,j]
                        dist = 1 - match/math.sqrt((size1 * size2))
                        rvars.append(dist)
                        avars.append( max( size1 ,  size2 ) -  match )

        return rvars,avars

def withoutvar(matrix, indexes1, indexes2):

        avars = []
        rvars = []
        for i in indexes1:

                size1 = max(1.0 , matrix[i,i])  

                #row = dm[i,...]

                for j in indexes2:

                        size2 = max(1.0, matrix[j,j])
                        match = matrix[i,j]

                        dist = 1 - match/math.sqrt((size1 * size2)) 
                        rvars.append(dist)
                        avars.append(  max( size1 ,  size2 ) -  match  )

        return rvars, avars


class node:

        def __init__(self, parent = None, name = "", distance =0.0, index = 0):

                self.children = []

                self.parent = parent

                self.genecount = ""

                self.name = name

                self.distance = distance

                self.mono = 1

                self.index = index

                self.annotation = 0

                if parent is not None:

                        parent.children.append(self)

        def __str__(self):

                if len(self.children) == 0:

                        return self.name+":"+str(self.distance)

                else:
                        return "("+",".join([str(x) for x in self.children])+")"+":"+str(self.distance)

        def push(self, name = "", distance =0.0, index = 0):

                newchild = node(self,name,distance,index)

                return newchild

        def build(self,text):

                if text == "":

                        return self

                elif text[-1] == ";":

                        text = text[:-1]

                current_node = self
                allnames = []

                index = 0
                current_name = ""
                for char in text:

                        if char in [" ", "\'"]:

                                continue

                        if char == "(":

                                current_node = current_node.push()

                        elif char == ")": 

                                name,distance = (current_node.name.split(":")+["0.0"])[:2]

                                name = name.split()[0].split(" ")[0].split("\t")[0]

                                if len(name)>0:

                                        current_node.name = name

                                        allnames.append(name)

                                        if len(current_node.children) == 0:
                                                current_node.index = index
                                                index += 1

                                else:

                                        allnames.append(current_node.name)

                                try:
                                        current_node.distance = float(distance.strip())

                                except:

                                        current_node.distance = 0.0

                                current_node = current_node.parent

                                if len(current_node.name) == 0 :

                                        current_node.name = "bh_"+str(len(allnames))

                                else:
                                        print(current_node.name)


                        elif char == "," or char ==";":

                                name,distance = (current_node.name.split(":")+["0.0"])[:2]

                                name = name.split()[0].split(" ")[0].split("\t")[0]

                                if len(name)>0:

                                        current_node.name = name

                                        allnames.append(name)

                                        if len(current_node.children) == 0:
                                                current_node.index = index
                                                index += 1

                                else:

                                        allnames.append(current_node.name)

                                try:
                                        current_node.distance = float(distance.strip())

                                except:

                                        current_node.distance = 0.0


                                if current_node.parent is not None:
                                        current_node = current_node.parent.push()

                        else:

                                current_node.name += char


                return self


        def allancestors(self):

                offset = 0.0
                ancestors=[]
                ancestor = self

                while ancestor.parent is not None:

                        ancestors.append(ancestor)

                        ancestor = ancestor.parent

                ancestors.append(ancestor)

                return ancestors


        def alloffsprings(self, offset = 0.0):

                alloffsprings = [(offset,self)]

                for child in self.children:

                        alloffsprings.extend(child.alloffsprings(offset = child.distance+offset))

                return alloffsprings

        def allleaves(self):
                return [x for x in self.alloffsprings() if len(x[1].children) == 0]
        def magnitude(self, exclude = None):

                if exclude is None or exclude is self:

                        alloffsprings = [x[0] for x in self.alloffsprings()]

                else:

                        exclude_offsprings = [x[1] for x in exclude.alloffsprings()]

                        alloffsprings = [x[0] for x in self.alloffsprings() if x[1] not in exclude_offsprings]


                return max(alloffsprings+[0.0])


        def findclades(self, dm , nametoclade):

                if len(self.children) == 0:
                        return []

                for child in self.children:
                        
                        child.findclades( dm , nametoclade)


                lnames = [x[1].name for x in self.children[0].allleaves()]
                rnames = [x[1].name for x in self.children[1].allleaves()]

                lclades = set([nametoclade[name] for name in lnames])
                rclades = set([nametoclade[name] for name in rnames])

                if len(lclades) != 1 or len(rclades) != 1 or lclades == rclades:
                        return

                lindexes = [x[1].index for x in self.children[0].allleaves()]
                rindexes = [x[1].index for x in self.children[1].allleaves()]


                lsize = sum([dm[i,i] for i in lindexes])/max(1,len(lindexes))
                rsize = sum([dm[i,i] for i in rindexes])/max(1,len(rindexes))


                without_rel,without_abs = withoutvar(dm, lindexes, rindexes)

                if len(lindexes) * len(rindexes)  < 10 or min( len(lindexes),  len(rindexes)) < 3:
                        
                        return 

                lwithin_rel,within_abs = withinvar(dm, lindexes)
                rwithin_rel,within_abs = withinvar(dm, rindexes)


                within_mean = ( sum(lwithin_rel) + sum(rwithin_rel) ) / max(len(lwithin_rel) + len(rwithin_rel), 1)
                without_mean = sum(without_rel) / len(without_rel)

                print(without_mean / within_mean if within_mean > 0 else 100)
                
                return 

        def splitbygenecount(self):

                if len(self.children) == 0:
                        return []

                clades = sum([child.splitbygenecount() for child in self.children], [])

                lindexes = [x[1].index for x in self.children[0].allleaves()]
                rindexes = [x[1].index for x in self.children[1].allleaves()]

                lcounts = set([x[1].genecount for x in self.children[0].allleaves()])
                rcounts = set([x[1].genecount for x in self.children[1].allleaves()])

                common = set.intersection(lcounts, rcounts)

                if len(lcounts) == 1 and len(common) == 0:

                        clades.append(self.children[0])

                if len(rcounts) == 1 and len(common) == 0:

                        clades.append(self.children[1])

                if len(common) == 1:
                        return [self]

                return clades


def distance_search(root, refs, samples,phyly_groups):

        if len(phyly_groups) == 0:

                return []


        refnodes = [root.search(ref) for ref in refs]

        results = []
        for group in phyly_groups:

                if group in refnodes:

                        continue

                if type(group) == type([]):

                        refnodes_within_group = [refnode[1] for node in group for refnode in node.get_offsprings(refs)]

                else:

                        refnodes_within_group = [x[1] for x in group.get_offsprings(refs)]

                if len(refnodes_within_group) < 1:

                        continue

                if type(group) == type([]):

                        result = [node for mono_group in group for node in mono_group.distinctive_nodes(refnodes_within_group) if len(node.allhaplotypes(samples))>0]

                else:

                        result = [node for node in group.distinctive_nodes(refnodes_within_group)]

                results.extend(result)


        dis_clusters = []
        for index,  mono in enumerate(results):

                mono_cluster = cluster("dis_"+str(index))

                mono_cluster.roots.append(mono)

                mono_cluster.members = [x[1] for x in mono.alloffsprings()]

                mono_cluster.label_members()

                dis_clusters.append(mono_cluster)


        return dis_clusters


def normtotdm(matrix):


        matrix = np.matrix(matrix)

        dm = np.ones((matrix.shape))
        vars = [matrix[i,i] for i in range(len(matrix))]
        for i in range(len(matrix)):

                var1 = max(1.0,vars[i])
                for j in range(i+1):

                        var2 = max(1.0,vars[j])

                        dm [i,j] = 1 - matrix[i,j]/math.sqrt((var1*var2))
                        dm [j,i] =  dm [i,j]

        return dm


def getgenecounts(genecountfile):

        nametoproteincoding = cl.defaultdict(list)
        with open(genecountfile, mode = 'r') as f:

                for line in f:

                        if len(line.strip()):

                                line = line.split('\t')

                                if  float(line[6]) > 90.0 and line[4] == "protein_coding":

                                        nametoproteincoding[line[0]].append(line[3].split('-')[0]) 

        nametoproteincoding = {name:";".join(sorted(genes)) for name,genes in nametoproteincoding.items()}


        return nametoproteincoding



def main(args):


        inputfile = args.input

        normfile = args.input + "_norm.txt"
        treefile = args.input + "_tree.ph"
        cladefile = args.input + "_types.txt"
        
        nametoclade = dict()
        with open(cladefile, mode ='r') as f:
                for line in f:
                        line = line.split('\t')
                        index,names = line[0],line[-1].strip()
                        for name in names.split(","):
                                nametoclade[name] = int(index)


        normmatrix = np.genfromtxt(normfile, delimiter=',')

        headers = dict()
        with open(args.input, mode = 'r') as f:
                for line in f:
                        if len(line) ==0 or line[0] != ">":
                                continue
                        lines = line.strip().split()
                        headers[lines[0][1:]] = lines[1]

        with open(treefile, mode = 'r') as f:

                tree_text = f.read()

        tree = node().build(tree_text)

        allclades = tree.findclades(normmatrix, nametoclade)

        


        return 0

def run():
        """
                Parse arguments and run
        """
        parser = argparse.ArgumentParser(description="program determine pangenome-allele types")
        parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
        parser.set_defaults(func=main)
        args = parser.parse_args()
        args.func(args)


if __name__ == "__main__":
        run()
