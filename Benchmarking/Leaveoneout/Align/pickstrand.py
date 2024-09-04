#!/usr/bin/env python3

allresults = "./all_results"

lastrow = [""]
with open(allresults, mode = 'r') as f:

        for line in f:

                line = line.split(",")

                if ".txt" in line[1]:
                        line = line[1:]

                if line[0].replace("_r.txt","") == lastrow[0].replace("_r.txt",""):

                        score1 = int(line[1])/max(1,int(line[2])+int(line[1]))
                        score2 = int(lastrow[1])/max(1,int(lastrow[2])+int(lastrow[1]))

                        if score1 < score2:

                                print(line[0].replace("_r.txt",""),int(line[1]), int(line[2]), score1)
                        else:
                                print(lastrow[0].replace("_r.txt",""),int(lastrow[1]), int(lastrow[2]),score2)

                lastrow = line
