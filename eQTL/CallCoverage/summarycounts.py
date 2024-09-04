#!/usr/bin/env python3

import os
import collections as cl

allfiles = ["results/"+x for x in os.listdir("results") if x.endswith("_count.txt")]
outputfile = "allcounts.txt"

samplenames = dict()
samplenamefile = "./Geuvadis_samples.txt"
with open(samplenamefile, mode = 'r') as f:
        for line in f:
                line = line.split()
                samplenames[line[1]] = line[0]

totalcounts = [0]*len(allfiles)
allgenes = cl.defaultdict(lambda:[0]*len(allfiles))
samples = [samplenames[thefile.split("/")[-1].split(".")[0]] for thefile in allfiles]
for index,thefile in enumerate(allfiles):

        with open(thefile, mode= 'r') as f:
                for line in f:
                        line = line.split()
                        allgenes[line[0]][index] = int(line[1])
                        totalcounts[index] += int(line[1])


with open(outputfile, mode = 'w') as f:

        f.write(",".join(["gene"]+samples)+"\n")

        for gene, counts in allgenes.items():

                counts = [str(round(x*1000000/total,5)) for x,total in zip(counts, totalcounts)]

                f.write(",".join([gene]+counts)+"\n")
