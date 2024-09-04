#!/usr/bin/env python3

import collections as cl

thefile = "summary.txt"

groups = cl.defaultdict(lambda: [0, 0, 0, 0])

with open(thefile, mode = 'r') as f:
        for line in f:
                line = line.strip().split()

                if int(line[1]) < 3:
                        continue

                groups [line[0]][1] += int(line[-1] == "1.0")
                groups [line[0]][0] += 1
                groups [line[0]][2] += int(line[1])
                if line[-1] == "1.0":
                        groups [line[0]][3] += int(line[1])

for groupname, count in groups.items():

        print(groupname, count[0], count[1], count[2],count[3])
