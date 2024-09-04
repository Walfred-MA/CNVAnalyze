#!/usr/bin/env python3

leave_allresults = "../allpairs/allpairs_strdpick"


genotype_result = "nonleave_genotype_match2.txt"

import collections as cl

"""
use_set = cl.defaultdict(lambda : 1.0)
with open(closest_allresults, mode = 'r') as f:
        for line in f:

                first, second = sorted(line.split(" ")[0].split("/")[-1].replace("_output.txt","").split("_pairwise_"))

                use_set[first] = min( use_set[first],float(line.split(" ")[-1]))
                use_set[second] = min(  use_set[second],float(line.split(" ")[-1]) )
"""

leave_results = cl.defaultdict(float)
with open(leave_allresults, mode = 'r') as f:
        for line in f:

                thepair = tuple(sorted(line.split(" ")[0].split("/")[-1].replace("_output.txt","").split("_pairwise_")))

                leave_results[thepair] = float(line.split(" ")[-1])

total = 0
nomismatch = 0  
lastrow = [""]
with open(genotype_result, mode = 'r') as f:

        for line in f:
                line = line.split(",")
                ori,genotype,closest  = line[1].split(";")[0],line[2].split(";")[0],line[3].split(";")[0]
                size = int(float(line[4]))


                thepair1 = tuple(sorted([ori,genotype]))
                thepair2 = tuple(sorted([ori,closest]))

                if ori==genotype or leave_results.get(thepair1,-1) == 0.0:
                        nomismatch+=1

                pair1result = 0.0 if thepair1[0] == thepair1[1] else leave_results.get(thepair1,-1)
                pair2result = 0.0 if thepair2[0] == thepair2[1] else leave_results.get(thepair2,-1)


                print(ori, pair1result, pair2result)

print(nomismatch)
