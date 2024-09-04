#!/usr/bin/env python3

import collections as cl
import argparse
import os
import re

def cigar_str(text):

        length = 0
        score = 0
        match = 0
        start = 0
        cigars = re.findall(r'\d+[A-Z]',text)

        for cigar in cigars:

                thesize, thetype = int(cigar[:-1]), cigar[-1]

                if thetype == 'S' or thetype == 'I':

                        length += thesize

                elif thetype == 'D':

                        start += thesize

                elif thetype == 'M':

                        start += thesize
                        match += thesize
                        length += thesize

        return match/length,start

class GenodeDB:

        def __init__(self):

                self.trans_exons = cl.defaultdict(list)
                self.exons_loci = cl.defaultdict(list)


        def add(self,chrom, start, end, name):

                self.exons_loci[chrom].append([start,end,name])

                return self

        def OverlapChr(self,  refs, queries):

                results = [0 for x in refs]

                length = len(refs) + len(queries)
                lgenes = len(refs)
                coordinates = [x for y in refs + queries for x in y[:2]] 

                coordinates_sortindex = sorted( range(len(coordinates)), key = lambda x: coordinates[x])


                current_refs = []
                current_queries = 0
                for index in coordinates_sortindex:

                        index2 = index // 2

                        if index % 2 == 0:
                                if index2 < lgenes:
                                        current_refs.append(index2)
                                        results[index2] += current_queries

                                else:
                                        current_queries += 1

                                        for ref in current_refs:
                                                results[ref] += 1
                        else:

                                if index2 < lgenes:
                                        current_refs.remove(index2)
                                else:
                                        current_queries -= 1


                return results

        def overlap(self, chrom, reads, genecounts):

                record = self.exons_loci[chrom]

                overlaps = self.OverlapChr(record, reads)

                for name, count in zip(record, overlaps):
                        if count:
                                genecounts[name[2]] += count

                return overlaps


def main(args):

        inputfile = args.input
        exonfile = args.exon

        exonDB = GenodeDB()
        with open(exonfile, mode = 'r') as f:

                for line in f:
                        if len(line.strip()) == 0:
                                continue

                        line = line.strip().split("\t")
                        name, exons, loci = line 

                        if len(loci) == 1:
                                continue
                        loci = loci.split(";")[:-1]

                        for locus in loci:
                                chrom, start, end = "_".join(locus.split("_")[:-2]), int(locus.split("_")[-2]), int(locus.split("_")[-1])

                                exonDB.add(chrom,start,end,name)

        genecounts = cl.defaultdict(int)
        lastchrom = ""
        coordinates = []
        index = 0
        with open(inputfile, mode = 'r') as f:

                for line in f:
                        if len(line.strip()) == 0:
                                continue

                        index += 1

                        if index % 1000000 == 0:
                                print(index)

                        line = line.strip().split()
                        if len(line) == 5:
                                chrom, start, flag, quality, cigar = line
                        else:
                                continue

                        score, end = cigar_str(cigar)

                        #overlaps = exonDB.overlap(chrom, int(start), int(start)+76)

                        if lastchrom and (chrom != lastchrom or index %1000000 == 0):

                                overlaps = exonDB.overlap(lastchrom, coordinates, genecounts )

                                coordinates.clear()

                        lastchrom = chrom

                        flag = bin(int(flag))[2:][::-1]

                        flag = str(flag)+'000000000000'[:12-len(flag)]

                        passfilter = int(flag[9] == '0' and flag[8]=='0' and flag[11]=='0' and flag[2]=='0' and flag[3]=='0')

                        if passfilter and score > 0.95:
                                coordinates.append([int(start), int(start) + end])

        outtable = []
        for gene, count in genecounts.items():

                outtable.append(gene+"\t"+str(count))


        with open (args.output, mode = 'w') as f:
                f.write("\n".join(outtable)+"\n")






def run():
        """
                Parse arguments and run
        """
        parser = argparse.ArgumentParser(description="program determine pangenome-allele types")
        parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
        parser.add_argument("-e", "--exon", help="path to input data file",dest="exon", type=str, required=True)
        parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
        parser.set_defaults(func=main)
        args = parser.parse_args()
        args.func(args)


if __name__ == "__main__":
        run()
