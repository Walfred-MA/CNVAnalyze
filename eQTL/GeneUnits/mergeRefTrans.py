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
                self.exons_totran = dict()
                self.transtype = cl.defaultdict(list)
                self.exonstogene = dict()
                self.transtogene = dict()
                self.selecttran = set()
                self.selectexon = set()
                self.selectgene = set()
                self.allexons = set()
                self.alltrans = set()

        def selectTran(self, apprisfile, gff3file):

                with open(apprisfile, mode = 'r') as f:

                        for line in f:
                                self.selecttran.add(line.split('\t')[2])

                gff3 = pd.read_csv(gff3file, header = None, sep = '\t')
                exons_texts = gff3 [(gff3[2] == 'exon' )].values.tolist()
                used_gene = set()
                transsize = cl.defaultdict(lambda: cl.defaultdict(int))
                for row in exons_texts:

                        exonstart,exonend,exonstrd = int(row[3]), int(row[4]),1 if row[6] == '+' else -1

                        text = row[-1]
                        text = text.split(';')

                        transcript_type  = [x for x in text if x.startswith("transcript_type=")][0][16:]
                        exon_id= [x for x in text if x.startswith("exon_id=")][0][8:]
                        if "." in exon_id:
                                exon_id = ".".join(exon_id.split(".")[:-1])
                        exon_num= int([x for x in text if x.startswith("exon_number=")][0][12:])
                        tran_id= [x for x in text if x.startswith("transcript_id=")][0][14:]
                        if "." in tran_id:
                                tran_id = ".".join(tran_id.split(".")[:-1])
                        gene_id= [x for x in text if x.startswith("gene_id=")][0][8:]

                        if "." in gene_id:
                                gene_id = ".".join(gene_id.split(".")[:-1])

                        if tran_id in self.selecttran:

                                used_gene.add(gene_id)

                        if "protein_coding" in row[-1] or "protein_id" in row[-1]:
                                transsize[gene_id][tran_id] += row[4] - row[3]


                for gene_id, record in transsize.items():

                        if gene_id in used_gene:
                                continue


                        largest_tran = sorted([(value, key) for key,value in record.items()],reverse = 1)[0][1]
                        self.selecttran.add(largest_tran)

                return self


        def readrecord(self, gff3file):

                gff3 = pd.read_csv(gff3file, header = None, sep = '\t')
                exons_texts = gff3 [(gff3[2] == 'exon' )].values.tolist()
                added_exons = set()

                for row in exons_texts:

                        exonstart,exonend,exonstrd = int(row[3]), int(row[4]),1 if row[6] == '+' else -1

                        text = row[-1]
                        text = text.split(';')

                        transcript_type  = [x for x in text if x.startswith("transcript_type=")][0][16:]
                        exon_id= [x for x in text if x.startswith("exon_id=")][0][8:]
                        if "." in exon_id:
                                exon_id = ".".join(exon_id.split(".")[:-1])
                        exon_num= int([x for x in text if x.startswith("exon_number=")][0][12:])
                        tran_id= [x for x in text if x.startswith("transcript_id=")][0][14:]
                        if "." in tran_id:
                                tran_id = ".".join(tran_id.split(".")[:-1])
                        gene_id= [x for x in text if x.startswith("gene_id=")][0][8:]

                        if "." in gene_id:
                                gene_id = ".".join(gene_id.split(".")[:-1])

                        self.allexons.add(exon_id)
                        if tran_id in self.selecttran:
                                self.selectexon.add(exon_id)


                        gene_name = [x for x in text if x.startswith("gene_name=")][0][10:]
                        tran_id = gene_name

                        self.alltrans.add(tran_id)

                        if exon_id not in added_exons:

                                self.exons_loci[row[0]].append([exonstart, exonend, exon_id])
                                self.transtogene[tran_id] = (gene_id, gene_name)
                                self.exonstogene[exon_id] = (gene_id, gene_name, row[4] - row[3], exonstrd * row[3], row[0])

                                self.exons_totran[exon_id] = (tran_id, len(self.trans_exons.get(tran_id, [])))
                                self.trans_exons[tran_id].append([exon_id, abs(row[4] - row[3])])
                                self.transtype[exon_id] = transcript_type
                                self.transtype[tran_id] = transcript_type

                                added_exons.add(exon_id)


                for chrom, record in self.exons_loci.items():

                        self.exons_loci[chrom] = sorted(record)


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

def overlap_scores(value):

        overlap_num = len([x for x in value if x[2] > 0.98* x[1] ])

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


                if (overlap_num1 or overlap_simi1 > 0.98) and  (overlap_num2 or overlap_simi2 > 0.98):

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

        gff3file = args.gff
        inputfile = args.input
        apprisfile = args.appris

        """
        refAlleles = cl.defaultdict(str)
        with open(table, mode = 'r') as f:
                for line in f:
                        line = line.split()
                        chrom = line[7].split(":")[0]

                        if chrom in mainchrom:
                                refAlleles[line[0]] = line[7]
        outpath = args.output
        """

        GffRecord = GenodeDB().selectTran(apprisfile,gff3file).readrecord(gff3file)

        transtogene = GffRecord.transtogene
        exons_totran = GffRecord.exons_totran
        trans_exons = GffRecord.trans_exons
        selectexon = GffRecord.selectexon

        table = pd.read_csv(inputfile, header = None, sep = '\t')
        table = table[(table[5].str.contains('ENSE') )]
        table = table.sort_values(by=[0,3], ascending=[1, 0])

        aligns = table[[0, 5, 6,2,3,4, 7,8,9,12]].values.tolist()

        aligns_bycontigs = cl.defaultdict(list)

        allgenenames = set()
        exonmatch_raw = cl.defaultdict(list)
        geneoverlaps = cl.defaultdict(int)
        for row in aligns:

                if '.' in row [1]:
                        row[1] = ".".join(row[1].split(".")[:-1])

                chunkname,exon,size,start,end,strand = row[:6]

                if exon not in GffRecord.exonstogene:
                        continue

                genename,genesize = GffRecord.exonstogene[exon][1],GffRecord.exonstogene[exon][2]

                allgenenames.add(genename)

                chrom,chunkstart,chunkend = "_".join(chunkname.split("_")[:-2]),int(chunkname.split("_")[-2]),int(chunkname.split("_")[-1])

                tran = GffRecord.exons_totran[exon]
                if strand == '-':
                        exonstart,exonend = chunkend - end, chunkend - start
                else:
                        exonstart,exonend = chunkstart + start, chunkstart + end

                overlaps = GffRecord.overlap(chrom,exonstart,exonend)

                keys = list(overlaps.keys())

                if row[-2] > 0.98*row[2] :

                        for exon1, overlap in overlaps.items():

                                osize = merge_segs(overlap)

                                genename1 = GffRecord.exonstogene[exon1][1]
                                if genename != genename1:
                                        geneoverlaps[(exon,exon1)] = max(geneoverlaps[(exon,exon1)],osize)

                if row[-2] > 0.9*row[2]:
                        for exon1, overlap in overlaps.items():

                                osize = merge_segs(overlap)
                                genename1 = GffRecord.exonstogene[exon1][1]
                                if genename != genename1 and osize > 0.9*genesize:

                                        exonmatch_raw[exon].append(exon1)

        #geneoverlaps = {key: merge_segs(value) for key,value in geneoverlaps.items()}

        transnonunique = cl.defaultdict(list)
        transoverlap = dict()
        exonmatch = cl.defaultdict(list)

        for pair, osize in geneoverlaps.items():

                exon1,exon2 = pair[0],pair[1]

                esize1,esize2 = GffRecord.exonstogene[exon1][2], GffRecord.exonstogene[exon2][2]

                tran1,tran2 = exons_totran[exon1],exons_totran[exon2]

                thekey = (tran1[0], tran2[0])

                if thekey not in transoverlap:

                        transoverlap[thekey] = [ [ x+[0] for x in trans_exons[tran1[0]] ],  [ x+[0] for x in trans_exons[tran2[0]] ]]

                transoverlap[thekey][0][tran1[1]][-1] = max(transoverlap[thekey][0][tran1[1]][-1], osize)
                transoverlap[thekey][1][tran2[1]][-1] = max(transoverlap[thekey][1][tran2[1]][-1], osize)

                if osize > 0.98 * esize1:

                        exonmatch[exon1].append(exon2)

        transoverlap_MANE = {key:[[x for x in value[0] if x[0] in selectexon],[x for x in value[1] if x[0] in selectexon]] for key,value in transoverlap.items()}

        allparas_index = mergepara(transoverlap_MANE)

        for genename in GffRecord.alltrans:
                if genename not in allparas_index:
                        allparas_index[genename] = max(list(allparas_index.values())+[-1]) + 1

        exonmatch_raw = {key:[x for x in value if allparas_index[GffRecord.exonstogene[x][1]] == allparas_index[GffRecord.exonstogene[key][1]]] for key,value in exonmatch_raw.items()}

        allexons_index = mergeexon(exonmatch_raw)

        for exon in GffRecord.allexons:
                if exon not in allexons_index:
                        allexons_index[exon] = max(list(allexons_index.values())+[-1]) + 1


        for tran in GffRecord.alltrans:

                transnonunique[tran] = [ [x[0],[]] for x in trans_exons[tran] ]

        allpairs = list(transoverlap.keys())
        for pair in allpairs:

                overlaps = transoverlap[pair]

                tran1,tran2 = pair[0], pair[1]

                if allparas_index[tran1] == allparas_index[tran2]:
                        continue

                for i,exon in enumerate(overlaps[0]):

                        genegroup = allparas_index[GffRecord.exonstogene[exon[0]][1]]

                        exonoverlap = [x for x in exonmatch[exon[0]] if allparas_index[GffRecord.exonstogene[x][1]] != genegroup ] 

                        transnonunique[tran1][i] = [exon[0], exonoverlap ]

                for i,exon in enumerate(overlaps[1]):

                        genegroup = allparas_index[GffRecord.exonstogene[exon[0]][1]]

                        exonoverlap = [x for x in exonmatch[exon[0]] if allparas_index[GffRecord.exonstogene[x][1]] != genegroup ] 

                        transnonunique[tran2][i] = [exon[0], exonoverlap] 

        transnonunique = {key:[x for x in value if x[0] in selectexon] for key,value in transnonunique.items()}

        exon_group_overlap = cl.defaultdict(list)
        for gene, exonoverlap in transnonunique.items():

                genegroup = allparas_index[gene]

                for exon in exonoverlap:

                        exon_index = allexons_index[exon[0]]

                        exon_group_overlap[exon_index].extend([x for x in exon[1] if allparas_index[GffRecord.exonstogene[x][1]] != genegroup])

        transnonunique_diffgroup = cl.defaultdict(list)

        for gene, exonoverlap in transnonunique.items():

                transnonunique_diffgroup[gene] = sorted([ [GffRecord.exonstogene[exon[0]][3], allexons_index[exon[0]], exon[0] ,exon_group_overlap[allexons_index[exon[0]]]] for exon in exonoverlap])


        genegroup = cl.defaultdict(list)
        for gene, index in allparas_index.items():

                genegroup[index].append(gene)

        outtable = []
        for index, group in genegroup.items():

                allunique_exons = []
                allloci = ""
                for gene in group:

                        unique_exons = [x[2] for x in transnonunique_diffgroup[gene] if len(x[-1]) == 0]
                        allunique_exons.extend(unique_exons)

                        loci = []
                        for exon in unique_exons:
                                loci.append([GffRecord.exonstogene[exon][4],abs(GffRecord.exonstogene[exon][3]),abs(GffRecord.exonstogene[exon][3])+ GffRecord.exonstogene[exon][2]])

                        allloci += ";".join(["_".join(map(str,x)) for x in loci])

                outtable.append([";".join(group)+";", ";".join(allunique_exons)+";",allloci +";"])


        pd.DataFrame.from_records(outtable).to_csv(args.output, index = False, header = None, sep = '\t', mode = 'w')

        return 



def run():
        """
                Parse arguments and run
        """
        parser = argparse.ArgumentParser(description="program determine psuedogene")
        parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
        parser.add_argument("-g", "--gff", help="path to input data file",dest="gff", type=str, required=True)
        parser.add_argument("-a", "--appris", help="path to input data file",dest="appris", type=str, required=True)
        parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)

        parser.set_defaults(func=main)
        args = parser.parse_args()
        args.func(args)


if __name__ == "__main__":
        run()
