randomAssignType = function(df)
{
        table = df[,-ncol(df)]
        colsums = colSums(table)
        probs <- colsums / sum(colsums)
        random_table <- matrix(0, nrow = nrow(table), ncol = ncol(table))
        colnum = ncol(table)
        # Generate random integers for each row
        for (i in 1:nrow(table)) {
                # Generate random integers based on row sum and column probabilities
                row_sum <- sum(table[i, ])  # Get row sum from the original table
                if (row_sum == 0)
                {
                        next
                }
                random_row <- sample.int(colnum, row_sum, replace = TRUE, prob = probs)
                # Assign the random integers to the corresponding row in the random table
                for (j in random_row)
                {
                        random_table[i, j] = random_table[i, j] + 1
                }
        }
        random_table = cbind(random_table, df[,ncol(df)])
        colnames(random_table) = colnames(df)

        return (data.frame(random_table))
}


randomAssign = function(df)
{
        table = df[,-ncol(df)]
        colsums = colSums(table)
        probs <- colsums / sum(colsums)
        random_table <- matrix(0, nrow = nrow(table), ncol = ncol(table))
        colnum = ncol(table)
        medianCN = median(rowSums(table))
        # Generate random integers for each row
        for (i in 1:nrow(table)) {
                # Generate random integers based on row sum and column probabilities
                row_sum <- medianCN  # Get row sum from the original table
                if (row_sum == 0)
                {
                        next
                }
                random_row <- sample.int(colnum, row_sum, replace = TRUE, prob = probs)
                # Assign the random integers to the corresponding row in the random table
                for (j in random_row)
                {
                        random_table[i, j] = random_table[i, j] + 1
                }
        }
        random_table = cbind(random_table, df[,ncol(df)])
        colnames(random_table) = colnames(df)

        return (data.frame(random_table))
}

randomSample = function(df)
{


        randomized_df <- df

        randomized_df[,ncol(df)] = df[,ncol(df)][sample(nrow(df), replace = FALSE)]

        return(randomized_df)

}
 # Define a function to run ANOVA model multiple times and compute average
trails <- function(model_formula, df, random, num_runs = 100) {
  # Initialize a vector to store results
  model_results <- c()
 
  # Run the model 'num_runs' times
  for (i in 1:num_runs) {
    # Generate randomized data frame
    randomized_df <- random(df)
    
    # Fit the ANOVA model
    model <- aov(model_formula, data = randomized_df)
    # Store the result (e.g., F-statistic, p-value, etc.)
    model_results <- c(model_results,sum(model$residuals^2))
  }
 
  # Calculate the average of the results
  average_result <- median(model_results)
  
  # Return the average result
  return(average_result)
}


args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]

table = read.table(input_file, sep = ',', header =T, row.names = 1)

df = data.frame(t(table))
response_var <- colnames(df)[ncol(df)]

fixed_effects <- c(colnames(df)[1:ncol(df)-1])
fixed_effects_str <- paste(fixed_effects, collapse = " + ")
formula_str <- paste(response_var, "~", fixed_effects_str, "")
model_formula <- as.formula(formula_str)


total_ss = var(df[,response_var]) * (nrow(df) -1 )
model <- aov(model_formula , data = df)
explained_variance = total_ss - sum(model$residuals^2)

model_randomtype_residuel = trails(model_formula, df, randomAssignType)
explained_ss_randomtype =  total_ss- model_randomtype_residuel

model_random_residuel = trails(model_formula, df, randomAssign)
explained_ss_random =  total_ss - model_random_residuel

model_random_residuel2 = trails(model_formula, df, randomSample)
explained_ss_random2 = total_ss - model_random_residuel2 



cat( input_file,nrow(df),ncol(df)-1  ,total_ss, explained_variance ,  explained_ss_random2, explained_ss_randomtype, explained_ss_random,  sep=',') 
