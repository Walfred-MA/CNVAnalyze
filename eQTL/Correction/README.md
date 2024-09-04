This is how to performed corrections on raw TPM counts from bamfiles.

PCA_fromPLINK2.sh: running plink2 command.

MergePCstoTPMs.py: the simple script to add first three PCs from plink2 to raw TPMs for correction.

correction_peer.R: the correction based on PEER for cross sample eQTL analysis, requiring first three PCs from reported chr1 genotypes.


correction_DESeq2.R: the correction based on DESeq2 for cross tissue eQTL analysis.
