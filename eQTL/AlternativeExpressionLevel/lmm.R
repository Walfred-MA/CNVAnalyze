library(car)

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

testallele = args[2]
#testallele = colnames(df)[testindex]

text <- args[3]
select_col <- strsplit(text, split = ",")[[1]]
#integer_vector <- as.integer(split_text)
select_col <- select_col[select_col %in% colnames(df) ]

if (length(select_col) == 0)
{
        quit()
}

df$nonalelle <- rowSums(df[, select_col])
table = df

express = table[,length(colnames(table))-1]
express = express - min(express)
#table = table[which(express < 2 * median(express )),]

response_var <- colnames(df)[length(colnames(table))-1]

#fixed_effects <- c(colnames(table)[3:length(colnames(table))-1])
fixed_effects <- c(colnames(df)[!colnames(df) %in% select_col])


fixed_effects <- fixed_effects[-length(fixed_effects)]
fixed_effects <- c(fixed_effects[-length(fixed_effects)], testallele ,"nonalelle")

# Create the formula string for fixed effects
fixed_effects_str <- paste(fixed_effects, collapse = " + ")

#SMN1 ~ exSMN1 + exSMN2 + exSMN0
# Complete formula string
formula_str <- paste(response_var, "~", fixed_effects_str, "+0")

# Convert to a formula object
model_formula <- as.formula(formula_str)

model <- lm(model_formula , data = table)
model = fix_alias(model)

formula_str2 <- paste(testallele, ' -  nonalelle = 0', collapse = "")
result = linearHypothesis(model, formula_str2 ,test='F')
summary(model)
as.data.frame(result)
