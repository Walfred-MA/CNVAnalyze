library(DESeq2)
library(BiocParallel)

args <- commandArgs(trailingOnly = TRUE)

expr = read.csv(args[1], header=TRUE, row.names = 1)

expr_matrix = as.matrix(expr)
expr_matrix <- apply(expr_matrix, 2, as.numeric)

colnames_split <- strsplit(colnames(expr_matrix), "\\.")
individual_ids <- sapply(colnames_split, `[`, 1)
tissue_types <- sapply(colnames_split, `[`, 2)

# Create the colData DataFrame
colData <- data.frame(
        individual_id = individual_ids,
        tissue_type = tissue_types,
        row.names = colnames(expr_matrix)
)

register(MulticoreParam(64)) 

dds <- DESeqDataSetFromMatrix(countData = expr_matrix, colData = colData, design = ~ tissue_type)

dds <- DESeq(dds)

normalized_counts <- counts(dds, normalized = TRUE)


geneindex <- rownames(normalized_counts)
genenames <- rownames(expr_matrix)

# Ensure that the order of gene names matches the order of the normalized counts
gene_names_aligned <- genenames[geneindex]

# Combine gene names and normalized counts into a single data frame
normalized_counts_df <- data.frame(GeneName = gene_names_aligned, normalized_counts)

# Write the normalized counts to file, including gene names and sample names
write.table(normalized_counts_df, file = "limma_normalized_counts.txt", sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
