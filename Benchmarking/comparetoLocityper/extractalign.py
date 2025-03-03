#!/usr/bin/env python3

import re
import os
import argparse


def alignmentcodi(seq,posi):
    if posi==[]:
        return []
    
    posi[0] +=1
    posi[1] -=1
   
     
    pieces =[a.span() for a in re.finditer(r'[A-Za-z]+',seq)]
    first=0
    if seq[0]=='-':
        first=1 
    s0=0
    index=0
    gaps=[]
    posifix=[a for a in posi]
    lastpos = 0
    totalgap = 0
    for piece in pieces:
        totalgap += piece[0] - lastpos
        s0=s0+piece[1] - piece[0]
        while index<len(posi) and posi[index]<s0:
            posifix[index]=posi[index]+totalgap
            index=index+1
        lastpos = piece[1]
    if index<len(posi):
        for i in range(index,len(posi)):
            posifix[index]=posi[index]+totalgap 
            
    posifix[0]-=1
    posifix[1]+=1
    
    return posifix

def findglobalalign(read, region1, region2):

    if read=='':
    
        return "",""

    alllines=[a.strip() for a in read.split('\n')]
    #ref: original genome, query:assemblies
    alllinesindex=[a for a in range(len(alllines)) if len(alllines[a])>0 and alllines[a][0]!='#' and alllines[a][0]!=':' and  re.match(r'\S+\s\S+',alllines[a])!=None]
    
    eachline=alllinesindex[::2]
    qlines=[a for a in eachline]
    rlines=[a+2 for a in eachline]
    alines=[a+1 for a in eachline]
    if qlines==[]:
        return "",""

    query=''.join([re.search(r'\s([A-Za-z]|-)+\s*',alllines[a]).group().strip() for a in qlines])
    ref=''.join([re.search(r'\s([A-Za-z]|-)+\s*',alllines[a]).group().strip() for a in rlines])
  
    
    acoordi1 = alignmentcodi(query,region1)
    
    qstart, qend = acoordi1
    
    acoordi2 = alignmentcodi(ref,region2)
    
    rstart, rend = acoordi2
   
    overlap_s = max(qstart, rstart) 
    overlap_e = min(qend, rend)

    return query[overlap_s:overlap_e].replace('-',''), ref[overlap_s:overlap_e].replace('-','')


def main(args):
    
    regions = args.region.split(",")
    region1 = [int(x) for x in regions[0].split("-")]
    region2 = [int(x) for x in regions[1].split("-")]
    
    with open(args.input, mode = 'r') as f:
        
        read = f.read()
        names = re.findall(r"# [12]:\s*(.+)", read )
        
    qalign, ralign = findglobalalign(read, region1, region2)
    
    with open(args.output, mode = 'w') as f:
        
        f.write(">{}\n{}\n>{}\n{}\n".format(names[0], qalign, names[1], ralign))
    
    


def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
    parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
    parser.add_argument("-r", "--region", help="path to input data file",dest="region", type=str, required=True)
    parser.add_argument("-o", "--output", help="path to input data file",dest="output", type=str, required=True)
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
