thefile = "../../fourthrun/PangenomeAlleles_type.tsv"

with open(thefile, mode ='r') as f:
        for line in f:

                line = line.split('\t')
                name, gene,rna,members = line[0], line[2],line[3],line[-1]


                if len(rna.strip()) ==0 or "(" in rna:
                        continue

                rna = { x.split(":")[0]:int(x.split(":")[-1]) for x in rna.split(";")}
                gene = gene.split(",")

                gene_new = [x for x in gene]
                for x in gene:
                        if x in rna:
                                gene_new.remove(x)
                                rna[x] -= 1


                themax = max([x for k,x in rna.items()] + [0])

                if themax > len(gene_new):
                        #print(line)
                        print(name,themax - len(gene_new), ",".join(gene) ,members.strip()+",")
