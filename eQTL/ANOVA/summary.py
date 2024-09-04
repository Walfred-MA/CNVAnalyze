#!/usr/bin/env python3

anovafile = "all_anova.txt"
snpanovafile = "all_snpanova.txt"
combanovafile = "all_combineanova.txt"

def mean_abs(thefile):

        with open(thefile, mode = 'r') as f:

                line = f.readline()
                allsums = [0] * len(line.strip().split(",")[1:])

                for line in f:

                        if "_group" not in line.split(",")[0]:
                                continue


                        line = list(map(float,line.strip().split(",")[1:]))

                        for i,value in enumerate(line):

                                allsums[i] += value


        mean_abs = []

        for i,x in enumerate(allsums):

                for j,y in enumerate(allsums[(i+1):]):

                        mean_abs.append(abs(x-y))

        return sum(mean_abs)/len(mean_abs), sum(allsums)/len(allsums)


records = dict()
with open(snpanovafile, mode = 'r') as f:

        for line in f:

                line = line.strip().split(",")
                name = line[0].split("_QTL_matrix.txt")[0].split("/")[-1]
                records[name] = line[1:]

records2  = dict()
with open(combanovafile, mode = 'r') as f:

        for line in f:

                line = line.strip().split(",")
                name = line[0].split("_Both_matrix.txt")[0].split("/")[-1]
                records2[name] = line[1:]





outtable = []
with open(anovafile, mode = 'r') as f:

        for line in f:

                line = line.strip().split(",")
                name = line[0].split("/")[-1]


                meanabs,mean_cn = mean_abs(line[0])

                record = records.get(name, ['0','0','0','0','0','0','0'])
                record2 = records2.get(name, ['0','0','0','0','0','0','0'])

                outtable.append(line + record+ [str(meanabs),str(mean_cn)]  + record2 )


with open(anovafile+"_merge.txt", mode = 'w') as f:

        f.write("\n".join([",".join(x) for x in outtable])+"\n")
