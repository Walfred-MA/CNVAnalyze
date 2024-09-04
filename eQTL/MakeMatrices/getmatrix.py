#!/usr/bin/env python3

import collections as cl

genotypefile = "../../types/PangenomeAlleles_typefix.tsv_Geuvadis_summary.txt"
PAannofile = "../All_PA.list"
#PAannofile ="../SMNSplice.list"
expressionfile = "../../peercorr/corrected_pheno_matrix.csv"
expressionfile = "../../bams/geuvadis/log/corrected_pheno_matrix.csv"
#expressionfile = "../../bams/geuvadis/allcounts.txt"
outfolder = "./matrix_log/"

PA_gene = cl.defaultdict(list)
with open(PAannofile, mode = 'r') as f:
        for line in f:
                if len(line.strip()) == 0:
                        continue

                line = line.strip().split()
                name = line[0]
                genes = line[-1].split(",")

                for gene in genes:
                        gene,size = gene.split(":")
                        if int(size) > 200:
                                PA_gene[name].append(gene)

type_gene = cl.defaultdict(list)
with open(genotypefile, mode = 'r') as f:

        header = f.readline()

        header = [x.split(".")[0] for x in header.split()[1:]]

        genotype_samples = {x:i for i,x in enumerate(header)}

        for line in f:
                if len(line.strip()) == 0:
                        continue

                line = line.strip().split()
                alleles = line[9]
                name = line[0]
                involvegene = set([x  for allele in alleles.split(",") for x in PA_gene[allele]])


                for gene in involvegene:
                        type_gene[gene].append([name]+line[10:])


with open(expressionfile, mode = 'r') as f:
        headers = f.readline().strip().split(",")[1:]
        headers = [x.split(".")[0] for x in headers]

        useheader = []
        headerindex = [-1] * len(headers)
        for i,header in enumerate(headers):

                if header in genotype_samples and headerindex[genotype_samples[header]] == -1:

                        headerindex[genotype_samples[header]] = i
                        useheader.append(header)

        for line in f:
                if len(line.strip()) == 0:
                        continue

                line = line.strip().split(",")
                gene = line[0][:-1]

                expressions = line[1:]
                expressions = [gene.replace(";","_")]+[expressions[i] for i in headerindex if i >= 0]

                if gene in type_gene:
                        type_gene[gene].append(expressions)

useheader = sorted(useheader, key = lambda x: genotype_samples[x])
for gene, data in type_gene.items():

        gene = gene.replace(";","_")
        with open(outfolder + gene +"_matrix.txt", mode = 'w') as f:

                line = ",".join([gene]+useheader)
                f.write(line+"\n")
                for line in data:

                        line = ",".join([str(x) for x in line])
                        f.write(line+"\n")
