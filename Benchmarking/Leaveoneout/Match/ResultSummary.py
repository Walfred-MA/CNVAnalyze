#!/usr/bin/env python3
import os
import collections as cl
import argparse
def getGenotypeResults(thefile, ResultsSummary):

        current_genegroup = ""
        current_sample = ""
        with open(thefile, mode = 'r') as f:

                for line in f:

                        if len(line) == 0:
                                continue

                        if line[0] == ">":

                                current_genegroup = line.split()[0][2:]
                                current_sample = line.split()[1].split('/')[-1].split('.')[0]   


                        elif line.startswith("result"):
                                line = line.split(":")[1]
                                ResultsSummary[current_genegroup][current_sample] = line.split(",")[:-1]



def main(args):

        inputfiles = [args.input+"/"+thefile for thefile in os.listdir(args.input) if thefile[-4:] == ".txt"]

        ResultsSummary = cl.defaultdict(lambda : dict())
        for thefile in inputfiles:

                getGenotypeResults(thefile, ResultsSummary)

        with open(args.output, mode = 'w') as f:

                for genegroup, genedata in ResultsSummary.items():

                        for sample, sampledata in genedata.items():

                                f.write("{}\t{}\t{}\n".format(genegroup,sample, ",".join(sampledata)))

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
