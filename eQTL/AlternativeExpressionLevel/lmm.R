library(car)

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


fix_alias <- function(model) {
        # Get the alias information
        alias_info <- alias(model)

        # Extract the names of the aliased coefficients
        aliased_coefs <- rownames(alias_info$Complete)

        if (length(aliased_coefs) == 0) {
                return(model)
        }

        # Get the original model formula
        model_formula <- formula(model)

        # Convert the formula to character for manipulation
        model_formula_char <- as.character(model_formula)

        # Extract the right-hand side of the formula
        rhs <- model_formula_char[3]

        # Split the rhs to get individual terms
        terms <- unlist(strsplit(rhs, " \\+ "))

        # Reverse the aliased_coefs to keep the last one
        aliased_coefs_reversed <- rev(aliased_coefs)

        # Remove redundant aliased coefficients while keeping the last one
        terms_to_keep <- setdiff(terms, aliased_coefs)

        # Reconstruct the formula without the redundant aliased coefficients
        new_rhs <- paste(terms_to_keep, collapse = " + ")
        new_formula <- as.formula(paste(model_formula_char[2], "~", new_rhs))

        # Fit the new model without redundant aliased coefficients
        new_model <- lm(new_formula, data = model$model)

        return(new_model)
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
#table = cleantable(df)
express = df[,length(colnames(df))]
#express = express - min(express)
#df = df[which(express < median(express ) + 2*sd(express)),]

if (length(rownames(df)) < 10)
{
        quit()
}

testalelle = args[2]
#testallele = colnames(df)[testindex]
response_var <- colnames(df)[length(colnames(df))]
text <- args[3]
select_col <- strsplit(text, split = ",")[[1]]
#integer_vector <- as.integer(split_text)
select_col <- select_col[which(select_col %in% colnames(df) & select_col != testalelle )]

df$nonalelle <- rowSums2(df[, select_col])

table = df

if (length(select_col) == 0 || (! testalelle %in% colnames(df)) || length(which(df$nonalelle > 0 ) ) <10 || length(which(df[,testalelle ] > 0 ) ) < 10  )
{
        quit()
}

df$nonalelle <- rowSums2(df[, select_col])
table = df
#express = express - min(express)

#fixed_effects <- c(colnames(table)[3:length(colnames(table))-1])
fixed_effects <- c(colnames(df)[!colnames(df) %in% select_col])


fixed_effects <- fixed_effects[-length(fixed_effects)]
fixed_effects <- c(fixed_effects[-length(fixed_effects)],"nonalelle")

# Create the formula string for fixed effects
fixed_effects_str <- paste(fixed_effects, collapse = " + ")
#SMN1 ~ exSMN1 + exSMN2 + exSMN0
# Complete formula string
formula_str <- paste(response_var, "~", fixed_effects_str, "+0")
# Convert to a formula object
model_formula <- as.formula(formula_str)

model <- lm(model_formula , data = table)
model = fix_alias(model)


# Extract coefficients
coefficients <- coef(model)


formula_str2 <- paste(testalelle, ' - nonalelle = 0', collapse = "")
result = linearHypothesis(model, formula_str2 ,test='Chisq')


testalelle_count = sum(table[,testalelle])
nonalelle_count = sum(table$nonalelle)

testalelle_express = coefficients[testalelle]
nonalelle_express = coefficients["nonalelle"]

cat(paste(nonalelle_count, nonalelle_express, sep = "\t"), "\n")
cat(paste(testalelle_count, testalelle_express, sep = "\t"), "\n")

as.data.frame(result)
