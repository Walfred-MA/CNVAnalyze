#!/usr/bin/env python3

import re
import pandas as pd
import os 
import argparse
import sys
import numpy as np
import multiprocessing as mul
import math
import collections as cl


def main(args):

        inputfile = args.input
        matrixfolder = args.matrix
        outputfolder = args.output

        mafiles = [x for x in os.listdir(matrixfolder) if x.endswith("_matrix.txt")]

        genetomfile = dict()
        for mfile in mafiles:

                genes = mfile.split(".")[0].split("_")
                for gene in genes:
                        genetomfile[gene] = matrixfolder + mfile




        names = []
        alleletogenes = dict()
        paralogs = cl.defaultdict(list)
        with open(inputfile, mode = 'r') as f:

                index = -1
                for line in f:

                        index += 1
                        line = line.split()
                        name, genes, locus = line[0], line[3],line[4]

                        genes = [x.split(":")[0].split("(")[0] for x in genes.split(";")] if genes != "NA" else []

                        alleletogenes[name] = genes 


                        for gene in genes:

                                paralogs[gene].append(index)
                        names.append(name)

        for refgene, group in paralogs.items():

                if len(group) < 2:
                        continue

                grouptext = ",".join([names[x] for x in group])

                for nameindex in group:

                        genes = alleletogenes[names[nameindex]]

                        for gene in [refgene]:
                                if gene in genetomfile:

                                        genefile = genetomfile[gene]
                                        prefix = genefile.split("/")[-1].split("_matrix.txt")[0]

                                        print("~/miniconda3/envs/eQTL/bin/Rscript lmm.R {} {} \"{}\" > {}/{}_result.txt ".format(genefile , names[nameindex], grouptext, outputfolder, prefix+"%"+names[nameindex]))   
                                        os.system("~/miniconda3/envs/eQTL/bin/Rscript lmm.R {} {} \"{}\" > {}/{}_result.txt ".format(genefile , names[nameindex], grouptext, outputfolder, prefix+"%"+names[nameindex]))



def run():
        """
                Parse arguments and run
        """
        parser = argparse.ArgumentParser(description="program determine pangenome-allele types")
        parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
        parser.add_argument("-m", "--matrix", help="path to input data file",dest="matrix", type=str, required=True)
        parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
        parser.set_defaults(func=main)
        args = parser.parse_args()
        args.func(args)


if __name__ == "__main__":
        run()
