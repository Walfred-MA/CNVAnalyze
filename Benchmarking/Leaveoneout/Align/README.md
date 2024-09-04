getrest.py: because too many pairs need to be aligned, so we always check for failed jobs and save successed jobs intead of finishing all of them once, this script is to find out the reserve jobs need to be done.

Snakemake_align_allpairs: snakemake pipeline to perform all alignments between each pair. 

pickstrand.py: the simple code to find the strands of each pair for comprison, in case they might in different strand.

comparetoleave.py: summary the results from alignments.

