#!/usr/bin/env python3

import collections as cl

thefile = "expr_summary.txt"
outfile = "expr_tissuetable.txt"

alldata = cl.defaultdict(lambda: cl.defaultdict(list))
with open(thefile, mode = 'r') as f:

        for line in f:

                if len(line.strip()) == 0:
                        continue

                line = line.strip().split()
                allele = line[2]

                gene = line[0]
                if allele != "nonallele":
                        alldata[gene][allele].append(lastline[0:5]+line[3:5]+[line[-1]])
                else:
                        lastline = line

outtable = []
for gene, genedata in alldata.items():

        for allele, alleledata in genedata.items():

                nonallelemax = sorted(alleledata, key = lambda x: float(x[5]) if x[5] != "NA" else -1000000, reverse = 1)[0]
                allelemax = sorted([x for x in alleledata if x[3] != nonallelemax[5]] , key = lambda x: float(x[3])  if x[3] != "NA" else -1000000, reverse = 1)[0]


                outtable.append([gene, allele, allelemax[1]] + allelemax[3:] + [nonallelemax[1]] +nonallelemax[3:])


with open(outfile, mode = 'w') as f:

        f.write("\n".join(["\t".join(x) for x in outtable]))
