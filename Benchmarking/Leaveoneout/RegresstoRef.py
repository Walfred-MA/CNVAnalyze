#!/usr/bin/env python3
import os
import collections as cl
import pandas as pd
import argparse
import math
import re

def normtodistm(norm, size):

        index = 0
        vars = []
        for i in range(size):
                vars.append(norm[index])
                index += (size - i)

        index = 0
        dm = [10000000000 for x in range(size**2)]
        for i in range(size):

                var1 = vars[i]

                for j in range(i, size):

                        var2 = vars[j]

                        themax = max(var1,var2)
                        #dm[i*size+j] = 1 - math.sqrt( (norm[index+j-i]**2)/(var1*var2+1) )
                        dm[j*size+i] = ( (themax - norm[index+j-i]), (norm[index+j-i])/max(themax, 1) )
                        dm[i*size+j] = dm[j*size+i]

                index += (size - i)


        return dm

def match(samplename, query_samples, ref_samples,  dm, samplegenes, vars):

        dm_size = len(samplegenes)

        mins = dict()
        for sample in query_samples:

                dises = dm[(dm_size*sample):(dm_size*sample+dm_size)] 

                themin_index = 0
                themin = (1000000000,0.0)
                for j, dis in enumerate(dises):

                        if j!=sample and dis < themin:

                                themin = dis
                                themin_index = j

                mins[sample] = (themin_index, themin)


        match_results = {}
        ref_lefted=set(ref_samples)

        query_samples = [x for x in query_samples if x not in match_results]
        ref_samples = [x for x in ref_samples ]

        ref_lefted=set(ref_samples)
        size = len(ref_samples)

        alldists = []
        for sample in query_samples:

                alldists.extend([min(dm[sample*dm_size + ref], dm[ref*dm_size + sample]) for ref in ref_samples])

        alldists_sortindex = sorted(range(len(alldists)), key = lambda x:alldists[x])
        for index in alldists_sortindex:

                query_index = query_samples[index//size]
                ref_index = ref_samples[index%size]

                if query_index in match_results:
                        continue

                theleast_value = alldists[index]

                match_results[query_index]=(ref_index,  theleast_value)

        for sample in query_samples:

                if sample not in match_results:

                        match_results[sample]= (-1,(0,0.0))

        leftindex = -1
        for sample in ref_lefted:
                match_results[leftindex] = (sample,(0,0.0))
                leftindex-=1

        results = []
        for query,ref in match_results.items():

                if query >= 0:
                        qname = samplegenes[query]
                        closestgene = samplegenes[mins[query][0]]
                        closestdis = mins[query][1]
                        qsize = vars[ query]
                        csize = vars[ mins[query][0]]

                else:
                        qname = samplename+"_NA"
                        closestgene = "NA"
                        closestdis = (0,0.0)
                        qsize =0 
                        csize = 0

                if ref[0] >= 0:
                        rname = samplegenes[ref[0]]
                        rsize = vars[ ref[0]]
                else:
                        rname = "NA"
                        rsize = 0

                results.append([samplename,qname, rname, closestgene, str(qsize) , str(ref[1][0]), str(ref[1][1]), str(rsize) , str(closestdis[0]),str(closestdis[1]), str(csize)])

        return results


def getClosest(regress_results,genenames,norm,vars):

        dist= normtodistm(norm, len(genenames))

        results = []
        for sample,regress in regress_results.items():

                regress_geneindex = [i for i,gene in enumerate(genenames)  if  re.search(r'chr[XY0-9]+_\d+',gene)]

                regress_set = set(regress)

                sample_geneindex = [[i]* len([x for x in regress if x in gene.split(";")]) for i,gene in enumerate(genenames)]
                sample_geneindex = [x for y in sample_geneindex for x in y ]


                results.extend(match(sample, sample_geneindex, regress_geneindex,  dist, genenames, vars))


        return results


def main(args):

        regress_results = cl.defaultdict(lambda:dict())
        with open(args.input, mode = 'r') as f:

                f.readline()

                for line in f:

                        if len(line) == 0:
                                continue
                        eles = line.split()+[""]

                        names = [x if eles[0] in x else eles[0] +"_"+x for x in eles[2].split(",")]

                        regress_results[eles[0]][eles[1]] = names




        lastmark = ''
        current_gene = ""
        genenames = []
        norm = []
        with open(args.output, mode = 'w') as w:
                with open(args.table, mode = 'r') as f:

                        for line in f:

                                if len(line) == 0:
                                        continue

                                if line[0] == "#":
                                        norm = []
                                        current_gene = line.split()[0][1:]
                                        genenames.clear()
                                elif line[0] == "%":
                                        kmersizes = list(map(float,line[1:].split()))   
                                elif line[0] == ">":
                                        names = ";".join([x.split()[0][1:] for x in line.split(";")])
                                        genenames.append(names)

                                elif line[0] == "$":

                                        norm += list(map(float,line[1:].split()))
                                elif lastmark == '$':

                                        results = getClosest(regress_results[current_gene],genenames,norm, kmersizes)
                                        for result in results:
                                                w.write(",".join(result)+"\n")
                                lastmark = line[0]


def run():
        """
                Parse arguments and run
        """
        parser = argparse.ArgumentParser(description="program determine psuedogene")
        parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
        parser.add_argument("-t", "--table", help="path to input data file", dest="table", type=str, required=True)
        parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str, required=True)
        parser.set_defaults(func=main)
        args = parser.parse_args()
        args.func(args)


if __name__ == "__main__":
        run()
