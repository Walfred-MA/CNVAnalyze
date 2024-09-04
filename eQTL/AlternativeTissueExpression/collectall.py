#!/usr/bin/env python3

import re
import os
import collections as cl

def extract_model_df(line):
        match = re.search(r'F-statistic: .+ on (\d+) and \d+ DF', line)
        if match:
                df_model = match.group(1)
                return df_model
        else:
                return "NA"

def readRoutput(file):

        estimated = []
        df = "NA"
        with open(file, mode = 'r') as f:

                lastline = ""
                for line in f:
                        if line.strip() == "---" or lastline.startswith("nonallele"):
                                estimated.append(lastline)
                        elif line.startswith("F-statistic:"):
                                df = extract_model_df(line)


                        lastline = line

        return estimated, df



results = "./results/"
outputfile = "expr_summary.txt"

allgenefiles = cl.defaultdict(list)
allfiles = [x for x in os.listdir(results) if x.endswith("_result.txt")]

for thefile in allfiles:

        gene = "_".join(thefile.split("_")[:-2])

        allgenefiles[gene].append(thefile)


Allresults = []
for gene, genefiles in allgenefiles.items():

        for genefile in genefiles:

                thetissue = genefile.split("_")[-2]

                Rresults,df = readRoutput(results+genefile)

                Allresults.extend([gene+"\t"+thetissue+"\t"+x.strip() + "\t" + df +"\n" for x in Rresults])

outtext = "".join([x for x in Allresults])+"\n"

with open(outputfile, mode = 'w') as f:

        f.write(outtext)
