#!/usr/bin/env python3

import pandas as pd
import os 
import argparse
import sys
import numpy as np
import multiprocessing as mul
import math
import collections as cl


def getdist(matrix, i ,j):

        return matrix[i,j]



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


class node:

        def __init__(self, parent = None, name = "", distance =0.0, index = 0):

                self.children = []

                self.parent = parent

                self.name = name

                self.distance = distance

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

        def alloffsprings(self, offset = 0.0):

                alloffsprings = [(offset,self)]

                for child in self.children:

                        alloffsprings.extend(child.alloffsprings(offset = child.distance+offset))

                return alloffsprings




def main(args):


        inputfile = args.input
        outputfile = args.output

        normfile = args.input + "_norm.txt"
        treefile = args.input + "_tree.ph"
        typefile = args.input + "_types.txt"

        nametogroup = dict()
        with open(typefile, mode = 'r') as f:
                lines = f.read().splitlines()
                for line in lines:
                        group = line.split(",")
                        for member in group:
                                nametogroup[member] = group



        normmatrix = np.genfromtxt(normfile, delimiter=',')

        matrix = normtotdm(normmatrix)

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

        allleaves = [x[1] for x in tree.alloffsprings() if len(x[1].children) == 0]

        nametoindex = {x.name:x.index for x in allleaves}
        allrefs = [(x.index,x.name) for x in allleaves if ("chr" in x.name and "v" not in x.name)]

        outputs = []
        for leave in allleaves:

                refdists = min([getdist(matrix, leave.index, ref[0]) for ref in allrefs])

                if len(nametogroup[leave.name])>2:
                        withintypedists =  sum([getdist(matrix, leave.index, nametoindex[member]) for member in nametogroup[leave.name] if member != leave.name]) /len(nametogroup[leave.name])
                else:
                        withintypedists = -1.0

                output =  [leave.name, refdists, withintypedists]

                thenode = leave 
                for level in range(8):

                        if thenode.parent is None:
                                output.append(-1.0)
                                continue

                        sibling = [x for x in thenode.parent.children if x != thenode ][0]
                        sibling_leaves = [x[1] for x in sibling.alloffsprings() if len(x[1].children) == 0]

                        sibling_leaves_dist = sum([getdist(matrix, leave.index, x.index) for x in sibling_leaves])/len(sibling_leaves)
                        output.append(sibling_leaves_dist)
                        thenode = thenode.parent

                outputs.append(output)


        #summary = [list(x) for x in zip(*outputs)]
        #summary = [sum([a for a in x if a >= 0 ])/len([a for a in x if a >= 0]) for x in summary[1:]]

        with open(args.output, mode = 'w') as f:

                f.write("\n".join(["\t".join([str(x) for x in output]) for output in outputs]) + "\n" )


        return 0

def run():
        """
                Parse arguments and run
        """
        parser = argparse.ArgumentParser(description="program determine psuedogene")
        parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
        parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str, required=True)
        parser.set_defaults(func=main)
        args = parser.parse_args()
        args.func(args)


if __name__ == "__main__":
        run()
