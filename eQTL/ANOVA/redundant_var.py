#!/usr/bin/env python3

import collections as cl


expressionfile = "../bams/geuvadis/log/redundant/corrected_pheno_matrix.csv"
outfile = "./redundant_var.txt"

with open(expressionfile, mode = 'r') as f:

        names = [x.split(".")[0] for x in f.readline().split(",")[1:]]
        nametoindex = dict()
        for name in names:
                nametoindex[name] = len(nametoindex)

        samesamples = cl.defaultdict(lambda: [[] for x in range(len(nametoindex)+1)] )
        for line in f:

                line = line.split(",")
                gene = line[0]
                values = [float(x) for x in line[1:]] 
                for name, value in zip(names,values):
                        samesamples[gene][nametoindex[name]].append(value)


with open(outfile, mode = "w") as f:
        for gene, pairs in samesamples.items():

                num = 0
                vars = 0.0
                for pair in pairs:
                        if len(pair) == 2:
                                vars += (pair[1]-pair[0])**2
                                num += 1

                if num:
                        f.write(gene+"\t"+str(vars/num)+"\n")
