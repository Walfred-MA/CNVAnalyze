rowSums2 <- function(x) {
        if (is.data.frame(x) || is.matrix(x)) {
                # For data frames or matrices with multiple columns
                row_sums <- rowSums(x)
        } else {
                # For single column vectors
                row_sums <- x
        }

        return(row_sums)
}




args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]

df = read.table(input_file, sep = ',', header =T)
df = as.data.frame(t(df))
# Set the first row as column names
colnames(df) <- df[1, ]
# Remove the first row
df <- df[-1, ]
df <- as.data.frame(lapply(df, function(x) as.numeric(as.character(x))))

colnames(df)[length(colnames(df))] = "express"


if (length(rownames(df)) < 10)
{
        quit()
}

multiple_model <- lm(express ~ . + 0, data = df)

data = df
data$aggreCN <- rowSums2(data[, colnames(data)[-length(colnames(data))]])



single_model <- lm(express ~ aggreCN + 0, data = data)

compare = anova(single_model, multiple_model)

f_value <- compare$"F"[2]  
p_value <- compare$"Pr(>F)"[2]  

cat(input_file, f_value, p_value ,'\n', sep = ",")
