#!/usr/bin/env python3

import pandas as pd
import collections as cl
import argparse
import pandas as pd

mainchrom = ['chr1', 'chr2','chr3','chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrUn']

class GenodeDB:

        def __init__(self):

                self.trans_exons = cl.defaultdict(list)
                self.exons_loci = cl.defaultdict(list)

        
        def add(self,chrom, start, end, name):
                
                self.exons_loci[chrom].append([start,end,name])
                
                return self

        def overlap(self,chrom, start, end):

                record = self.exons_loci[chrom]

                size = abs(end - start)
                overlaps = cl.defaultdict(list)

                for x in record:

                        second_start = max(start, x[0])
                        first_end = min(end, x[1])
                        if second_start - first_end < -20:

                                overlaps[x[2]].append([second_start,first_end])

                return overlaps

def merge_segs(segments):

        segments = [x for y in segments for x in y[:2]]

        segments_sortindex = sorted(range(len(segments)), key = lambda x: segments[x])

        start = 0
        merged = 0
        depth = 0
        for i, x in enumerate(segments_sortindex):


                if x %2:
                        depth -= 1
                        if depth == 0:
                                merged += segments[x] - start
                else:
                        depth += 1
                        if depth == 1:
                                start = segments[x]

        return merged

def getgroupoverlap(overlap, exontogroup):
        
        groupoverlaps = cl.defaultdict(list)
        
        for exon, overlap in  overlap.items():
                
                exongroup = exontogroup[exon]
                groupoverlaps[exongroup].extend(overlap)
                        
        groupoverlap = {exongroup:merge_segs(overlaps) for exongroup, overlaps in groupoverlaps.items()}
        
        return groupoverlap
                

def overlap_scores(value):

        overlap_num = len([x for x in value if x[2] > 0.97* x[1] ])

        overlap_size1 = max(1,sum([x[1] for x in value]))
        overlap_size2 = sum([x[2] for x in value])

        return int(overlap_num == len(value) and overlap_num > 0), overlap_size2/overlap_size1


def mergeexon(transoverlap):


        allexons =  sorted(list(set(list(transoverlap.keys()) + sum(list(transoverlap.values()),[]))))
        allexons = {exon:i for i,exon in enumerate(allexons)}
        allexons_index = list(range(len(allexons)))

        for exon, exonoverlaps in transoverlap.items():

                exonoverlaps = exonoverlaps + [exon]

                themin = min([allexons_index[allexons[exon]] for exon in exonoverlaps])

                for exon in exonoverlaps:

                        allexons_index[allexons[exon]] = themin

        for index in list(range(len(allexons))):

                newindex = allexons_index[index] 

                while newindex != index:

                        index = newindex
                        newindex = allexons_index[index] 

                allexons_index[index] = newindex

        allexons = {exon:allexons_index[i] for exon,i in allexons.items()}

        return allexons


def mergepara(transoverlap):

        allparas = sorted(list(set([x for y in transoverlap.keys() for x in y])))
        allparas = {para:i for i,para in enumerate(allparas)}
        allparas_index = list(range(len(allparas)))

        allmatches = set()
        for key, values in transoverlap.items():

                key1,key2 = key
                exon1,exon2 = values

                overlap_num1, overlap_simi1 = overlap_scores(exon1)
                overlap_num2, overlap_simi2 = overlap_scores(exon2)


                if (overlap_num1 or overlap_simi1 > 0.97) and  (overlap_num2 or overlap_simi2 > 0.97):

                        index1 = allparas[key1]
                        index2 = allparas[key2]

                        allparas_index[index1] = min(allparas_index[index1],allparas_index[index2])
                        allparas_index[index2] = min(allparas_index[index1],allparas_index[index2])


        for index in list(range(len(allparas))):

                newindex = allparas_index[index] 

                while newindex != index:

                        index = newindex
                        newindex = allparas_index[index] 

                allparas_index[index] = newindex


        allparas = {para:allparas_index[i] for para,i in allparas.items()}

        return allparas




def main(args):

        inputfile = args.input
        annofile = args.anno
        allelefile = args.title
                        
        exontogroup = dict()
        with open(annofile, mode = 'r') as f:
                for line in f:
                        gname, exons, loci = line.strip().split("\t")
                        
                        for exon in exons.split(";")[:-1]:
                                exontogroup[exon] = gname[:-1]
       
        try: 
                table = pd.read_csv(inputfile, header = None, sep = '\t')
        except:
                try:
                        os.system("touch " + args.output)
                except:
                        pass

                return 

        table = table[(table[5].str.contains('ENSE') )]
        table = table.sort_values(by=[0,3], ascending=[1, 0])
        aligns = table[[0, 5, 6,2,3,4, 7,8,9,12]].values.tolist()
        
        exonloci = GenodeDB()
        for row in aligns:
                
                if '.' in row [1]:
                        row[1] = ".".join(row[1].split(".")[:-1])
                        
                chunkname,exon,size,start,end,strand = row[:6]
                
                if exon not in exontogroup or row[-2] <= 0.98*row[2]:
                        continue
                
                chrom,chunkstart,chunkend = "_".join(chunkname.split("_")[:-2]),int(chunkname.split("_")[-2]),int(chunkname.split("_")[-1])
                
                exonstart,exonend = chunkstart + start, chunkstart + end
        
                exonloci.add(chrom, exonstart,exonend , exon)
        

        output = []
        with open(allelefile, mode = 'r') as f:
                for line in f:
                        
                        if len(line.strip()) == 0:
                                continue
                        
                        line = line.split()
                        
                        name, locus = line[0][1:], line[1]
                        
                        chrom, coordi = locus.split(":")[0], locus.split(":")[1][:-1].split("-")
                        
                        start,end = int(coordi[0]), int(coordi[1])
                        overlap = exonloci.overlap(chrom,start,end)

                        groupoverlap = getgroupoverlap(overlap, exontogroup)
                        
                        groupoverlap = ",".join([(k+":"+str(v)) for k,v in groupoverlap.items()])
                        
                        if len(groupoverlap):
                                output.append("\t".join([name, locus, groupoverlap]))

        with open(args.output, mode = 'a') as f:
                
                f.write("\n".join(output)+"\n")


def run():
        """
                Parse arguments and run
        """
        parser = argparse.ArgumentParser(description="program determine psuedogene")
        parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
        parser.add_argument("-a", "--annotate", help="path to input data file",dest="anno", type=str, required=True)
        parser.add_argument("-t", "--title", help="path to output file", dest="title",type=str, required=True)
        parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)

        parser.set_defaults(func=main)
        args = parser.parse_args()
        args.func(args)


if __name__ == "__main__":
        run()
