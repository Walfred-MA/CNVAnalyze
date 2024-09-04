#!/usr/bin/env python3


thepool = "./allpairs_strdpick"

finishpairs = set()
with open(thepool, mode = 'r') as f:
        for line in f:
                line = line.split()[0].split("/")[-1]
                thepair = line.split("_pairwise_")
                thepair[1] = thepair[1].split("_output.txt")[0]
                finishpairs.add(tuple(sorted(thepair)))

matchfile = "../leaves/leave_genotype_match2.txt"
matchfile2 = "../nonleaves/nonleave_genotype_match2.txt"
matchfile3 = "../refer/ref_genotype_match.txt"
allpairs = set()

with open(matchfile, mode = 'r') as f:
        for line in f:
                line = line.split(",")
                if "NA" in line:
                        continue

                allele, genotye, closest = line[1].split(";")[0],line[2].split(";")[0],line[3].split(";")[0]
                pair1 = tuple(sorted([allele, genotye]))
                pair2 = tuple(sorted([allele, closest]))


                allpairs.add(pair1)
                allpairs.add(pair2)


with open(matchfile2, mode = 'r') as f:
        for line in f:
                line = line.split(",")
                if "NA" in line:
                        continue

                allele, genotye, closest = line[1].split(";")[0],line[2].split(";")[0],line[3].split(";")[0]
                pair1 = tuple(sorted([allele, genotye]))
                pair2 = tuple(sorted([allele, closest]))


                allpairs.add(pair1)
                #allpairs.add(pair2)

with open(matchfile3, mode = 'r') as f:
        for line in f:
                line = line.split(",")
                if "NA" in line:
                        continue

                allele, genotye, closest = line[1].split(";")[0],line[2].split(";")[0],line[3].split(";")[0]
                pair1 = tuple(sorted([allele, genotye]))
                pair2 = tuple(sorted([allele, closest]))


                allpairs.add(pair1)
                #allpairs.add(pair2)

rest = allpairs - finishpairs

outfile = "./restpairs.txt"


with open(outfile , mode = 'w') as f:
        for pair in rest:

                if pair[0] != pair[1]:
                        f.write(",".join(list(pair))+"\n")
