#!/usr/bin/env python3

import os
import collections as cl
import argparse
import re


class trimer:
    
    def __init__(self, file):
        
        self.regions = cl.defaultdict(list)
        with open(file, mode = 'r') as f:
            
            for line in f:
                
                line = line.split()
                self.regions[line[0]].append([int(line[1]),int(line[2])])
        

    def trim(self, read, chr, start, end, strd):
       
 
        regions = self.regions[chr]
        size = end - start
        
        overlaps = [region for region in regions if max(end, region[1]) - min(start, region[0]) - size - (region[1] - region[0]) < 0]
       
        overlaps = sorted(overlaps, key = lambda x: x[0]-x[1])
        if len(overlaps):
 
            region = overlaps[0]

            overlap = [max(start, region[0]),min(end, region[1])] 
           

            if strd == "+":
                relative_s = overlap[0] - start
                relative_e = overlap[1] - start
            else:
                relative_s  = end - overlap[1]
                relative_e  = end - overlap[0]
           
 
            return read[relative_s:relative_e],relative_s,relative_e
        
        else:
            
            return "",0,0
        
    
def main(args):
   
    thetrimer = trimer(args.region)
    
    with open(args.input, mode = 'r') as f:
        
        reads = [x.splitlines() for x in f.read().split(">")[1:]]
        
        reads = [[read[0],"".join(read[1:])] for read in reads]
        
    
    with open(args.output, mode = 'w') as f:
        
        for read in reads:
         
            name, region = read[0].split()[:2]
            strd = region[-1]
            chr, region = region[:-1].split(":")
            region = region.split("-")
            start,end = int(region[0]),int(region[1])
 
            trim,relative_s,relative_e = thetrimer.trim(read[1],chr,start, end, strd)
           
            if len(trim ):
            
                f.write(">{}\t{}\t{}\n{}\n".format(name, relative_s,relative_e , trim))
            
            
def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program determine psuedogene")
    parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
    parser.add_argument("-o", "--output", help="path to input data file", dest="output", type=str, required=True)
    parser.add_argument("-r", "--region", help="path to input data file", dest="region", type=str, required=True)
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
