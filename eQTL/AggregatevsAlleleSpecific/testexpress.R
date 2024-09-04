args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]

df = read.table(input_file, sep = ',', header =T)
df = as.data.frame(t(df))
# Set the first row as column names
colnames(df) <- df[1, ]
# Remove the first row
df <- df[-1, ]
df <- as.data.frame(lapply(df, function(x) as.numeric(as.character(x))))
#table = cleantable(df)
express = df[,length(colnames(df))]
#express = express - min(express)

if (length(rownames(df)) < 10)
{
        quit()
}

aggreCN = rowSums(df) - df[,length(colnames(df))]
express = df[,length(colnames(df))]

spearman_test <- cor.test(aggreCN, express, method = "pearson")

# Extract the correlation coefficient and p-value
spearman_correlation <- spearman_test$estimate
spearman_p_value <- spearman_test$p.value

cat(input_file, spearman_correlation, spearman_p_value,'\n', sep = ",")
