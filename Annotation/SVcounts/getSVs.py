#!/usr/bin/env python3

import collections as cl

thefile = "../../fourthrun/PangenomeAlleles_annotation.tsv"

SVbyclass = cl.defaultdict(list)

with open(thefile, mode = 'r') as f:

        for line in f:

                line = line.split("\t")

                name = line[0]  
                ifintron = int(line[6] != '1')
                groupname = name.split('_')[0]

                if  line[6] in ["Intron","Decoy"]:
                        continue

                SVtype = "NA"

                name, gene,locus, rna,theclass, tag, thesv = line[0], line[2],line[3], line[4], line[5], line[13], line[14]

                if len(rna.strip()) and "(" not in rna and theclass!= "Ref":

                        rnas = cl.defaultdict(int)
                        for x in rna.split(";")[:-1]:

                                rnas[x.split(":")[1]] += int(float(x.split(":")[-1])*0.01+0.1)

                        gene_ = gene.split(",")
                        gene_new = [x for x in gene_]
                        for x in gene_:
                                if x in rnas:
                                        gene_new.remove(x)
                                        rnas[x] -= 1


                        themax = max([x for k,x in rnas.items()] + [0])

                        if themax > len(gene_new):

                                SVtype = "FullgeneDup"

                                if themax - len(gene_new) > 2:
                                        SVtype+="RunAway"

                                print(name,gene,groupname+"_"+locus, theclass.replace("Mutate","Alt"), SVtype)
                                continue



                if thesv != "NA":

                        SVs = [x.split("_") for x in thesv.split(";")]
                        SVs = sorted([(int(x[-1]),x[0]) for x in SVs],reverse = 1)[0]
                        if SVs[0] > 300:
                                SVtype = SVs[1]

                        if "Contraction" in SVtype or "Dele" in SVtype:
                                SVtype = "Deletion"

                        elif "Insert" in SVtype:
                                SVtype = "Insertion"

                        if SVtype == "Deletion" and "Exon" in tag:
                                SVtype += "Exon"

                ifconvert = "Not"
                if line[8] != line[11] and line[11] != "NA":
                        ifconvert = "Convertion"


                print(name,gene,groupname+"_"+locus,theclass.replace("Mutate","Alt"), SVtype, ifconvert)
