#!/usr/bin/env python3

typefile = "/panfs/qcb-panasas/wangfeim/Build/fourthrun/PangenomeAlleles_typefix.tsv"

outputfile  = "allpvalues.txt"

inputfile = "summary.txt"

typename = dict()
with open(typefile , mode = 'r') as f:
        for line in f:
                line = line.strip().split()
                typename["_".join(line[0].split("_")[:-1])] = [line[1], line[4].split(":")[0]]


outtext = []
with open(inputfile , mode = 'r') as f:
        for line in f:
                line = line.strip().split()

                outtext.append("\t".join([line[0],typename["_".join(line[1].split("_")[:-1])][0],typename["_".join(line[1].split("_")[:-1])][1]]  + line[1:])+"\n")


with open(outputfile, mode = 'w') as f:
        f.write("\n".join(outtext))
