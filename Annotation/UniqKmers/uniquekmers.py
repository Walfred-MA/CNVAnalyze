#!/usr/bin/env python3


import re
import os
import collections as cl
import argparse

def intdecode(string):

        number = 0

        for i,c in enumerate(string):

                code = ord(c) - ord('0')
                number *= 64
                number += code 

        return number


def main(args):

        matrixfile = args.matrix

        typefile = args.type


        with open(typefile, mode = 'r') as f:
                groups = [line.strip().split('\t')[-1].split(',') for line in f.readlines() if len(line.strip())]


        groupindex= dict()
        for i,group in enumerate(groups):

                groupindex.update({x:i for x in group})


        group_byindex = [set([ ]) for x in groups]
        results =  []
        kmerindex = -1
        kmernum = 0
        nameindex = 0

        group_counts = [0 for x in groups]
        ifunique = [0 for x in groups]
        with open(matrixfile, mode = 'r') as f:

                for line in f:

                        if len(line.strip()) == 0:
                                continue

                        if line[0] == '>':

                                names = [x.split()[0][1:] for x in line.split(';')]
                                for name in names:

                                        group_byindex[groupindex[name]].add(nameindex)
                                        group_counts[groupindex[name]]+=1
                                nameindex += 1


                        if line[0] == '$':
                                break

                allindex = set(range(nameindex))

                for line in f:


                        if len(line.strip()) == 0 or line[0] != '&':
                                continue

                        kmerindex += 1
                        kmernum += 1

                        if line[1] not in  ['-', '+'] :
                                continue

                        kmersamples =  set([ intdecode(x) for x in line.split()[-1].split(',')[:-1]])

                        if line[1] == '-':

                                kmersamples = allindex - kmersamples

                        for i,group in enumerate(group_byindex):

                                if set(kmersamples) == group:

                                        ifunique[i] = 1

                                union = kmersamples.union(group)
                                common = kmersamples.intersection(group)

                                ratio = len(common)/ len(union)

                                ifunique[i] = max(ratio, ifunique[i])


        [print(str(s)+'\t'+str(x)) for s,x in zip( group_counts  ,ifunique)]






def run():
        """
                Parse arguments and run
        """
        parser = argparse.ArgumentParser(description="program compiles sequences to kmer matrix")
        parser.add_argument("-m", "--matrix", help="path to input data file", dest="matrix", type=str, required = True)
        parser.add_argument("-t", "--type", help="path to input data file", dest="type", type=str, required = True)
        parser.set_defaults(func=main)
        args = parser.parse_args()
        args.func(args)


if __name__ == "__main__":
        run()
