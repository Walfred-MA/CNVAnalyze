./plink2 --vcf input.vcf.gz --make-bed --out input

./plink2 --bfile input --maf 0.05 --indep-pairwise 50 5 0.2 --out pruned_data

./plink2 --bfile input --extract pruned_data.prune.in --make-bed --out pruned_input

./plink2 --bfile pruned_input --pca 3 --out pca_results
