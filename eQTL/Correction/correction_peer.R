#!/usr/bin/env python3

library(peer)

args <- commandArgs(trailingOnly = TRUE)

expr = read.csv(args[1], header=TRUE, row.names = 1)
PCs = data.frame(read.csv(args[2], header=TRUE, row.names = 1, sep = '\t'))

samplenames = colnames(expr)
#samplenames<- as.vector(sapply(samplenames, function(x) paste(strsplit(x, "\\.")[[1]][1], collapse = ".")))

samplenames <- as.vector(sapply(samplenames, function(x) {
  parts <- strsplit(x, "\\.")[[1]]
  if (length(parts) >= 2) {
    paste(parts[1:2], collapse = ".")
  } else {
    parts[1]
  }
}))


valid_columns <- samplenames %in% colnames(PCs)
expr <- expr[,valid_columns]
covariates <- PCs[c("PC1", "PC2", "PC3"),samplenames[valid_columns]]


expr_matrix = as.matrix(expr)
expr_matrix <- apply(expr_matrix, 2, as.numeric)
expr_matrix = log2(expr_matrix + 1)
expr_matrix = t(expr_matrix)
model = PEER()
PEER_setPhenoMean(model,expr_matrix)
PEER_setNk(model,3)
dim(PEER_getPhenoMean(model))
# Retrieve the phenotype matrix
PEER_update(model)

# Retrieve the corrected phenotype matrix (residuals)
corrected_matrix <- PEER_getResiduals(model)

corrected_matrix = t(corrected_matrix)

corrected_matrix = 2^corrected_matrix
# Set the column names and row names to match those of expr
colnames(corrected_matrix) <- colnames(expr)
rownames(corrected_matrix) <- rownames(expr)


# Combine row names as the first column for writing to CSV
output_matrix <- cbind(RowIndex = rownames(corrected_matrix), corrected_matrix)

# Write the matrix to a CSV file without quoting numeric values
write.csv(output_matrix, file = "corrected_pheno_matrix.csv", row.names = FALSE, quote = FALSE)

# Confirm the dimensions of the corrected phenotype matrix
factors <- PEER_getX(model)
rownames(factors) <- colnames(expr)
write.csv(factors, file = "corrected_pheno_factors.csv", row.names = TRUE, quote = FALSE)

weights <- PEER_getW(model)
#colnames(weights) <- colnames(expr)
write.csv(weights, file = "corrected_pheno_weights.csv", row.names = TRUE, quote = FALSE)

