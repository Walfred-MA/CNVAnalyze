#!/usr/bin/env python3

import os
import pandas as pd

file = "allrst.txt"
outfile = "All.geuv_coordinates.txt"
liftoverfolder = "/project/mchaisso_100/cmb-16/walfred/Tools/"

index = 0

run = 0
if run:
        with open(file, mode = 'r') as r:

                with open(outfile, mode = 'w') as w:

                        for line in r:

                                if len(line.strip())==0:
                                        continue

                                line = line.split()

                                if not line[1].startswith("chr"):
                                        coordis = [x[1:-1] for x in line[4:]]
                                else:
                                        coordis = [line[1]]

                                for coordi in coordis:
                                        if coordi.startswith("chr"):

                                                chrom, posi = coordi.split(".")
                                                index += 1
                                                w.write("{}\t{}\t{}\t{}\n".format(chrom, posi ,posi, index))

#os.system("{}/liftOver {} {}/hg19ToHg38.over.chain {}_tohg38.txt {}_tohg38.txt_unmap.txt ".format(liftoverfolder, outfile, liftoverfolder, outfile, outfile))

os.system("{}/liftOver {} {}/hg38ToHg19.over.chain {}_tohg19.txt {}_tohg19.txt_unmap.txt ".format(liftoverfolder, outfile, liftoverfolder, outfile, outfile))

# Load original and lifted coordinates
original = pd.read_csv(outfile, sep='\t', header=None, names=['chr_orig', 'start_orig', 'end_orig', 'name'])
lifted = pd.read_csv(outfile+"_tohg19.txt", sep='\t', header=None, names=['chr_lifted', 'start_lifted', 'end_lifted', 'name'])

# Merge on the name column
combined = pd.merge(original, lifted, on='name')

# Reorder columns if necessary
combined = combined[['chr_orig', 'start_orig', 'end_orig', 'name', 'chr_lifted', 'start_lifted', 'end_lifted']]

# Save to a new file
combined.to_csv(outfile+"_tohg19.txt", sep='\t', header=False, index=False)
