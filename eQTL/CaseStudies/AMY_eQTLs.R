
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

table <- read.csv("~/Documents/Marklab/eQTL/AMY2B_matrix.txt", header=FALSE, skip = 1, row.names=1)
table <- t(table)
table = data.frame(table)
table$express = table[,"AMY2B"]

model = lm(express ~ AMY_group1_0_0 + AMY_group1_1_0 + AMY_group1_2_0 + AMY_group1_3_0 + AMY_group1_4_1 + AMY_group1_34_11, data = df)

AMY2B_ = "AMY_group1_0_0,AMY_group1_1_0,AMY_group1_2_0,AMY_group1_3_0"
Trans  = "AMY_group1_4_1,AMY_group1_34_11"

AMY2B_ = strsplit(AMY2B_, split = ",")[[1]]
Trans  = strsplit(Trans, split = ",")[[1]]

table$allAMY2B = rowSums2(table[, AMY2B_])
table$allTrans  = rowSums2(table[, Trans ])

model_formula <- as.formula("express ~ allAMY2B  + allTrans")

model <- lm(model_formula , data = table)

result = linearHypothesis(model, "allAMY2B - allTrans" ,test='Chisq')
result

summary(model)

model_formula <- as.formula("express ~ allAMY2B  + allTrans + 0")

model <- lm(model_formula , data = table)

result = linearHypothesis(model, "allAMY2B - allTrans" ,test='Chisq')
result

summary(model)


x = unlist(table$allAMY2B)
y = unlist(table$AMY_group1_4_1)
z = unlist(table$allTrans)

fourth_dim = unlist(table$express)
datafull <- data.frame(x, y, z, fourth_dim)

express_null = datafull [which(datafull $x==2 & datafull$z == 0),"fourth_dim"]
express_amy2b = datafull [which(datafull$x > 2 & datafull$z == 0),"fourth_dim"]
express_tran1 = datafull [which(datafull$y == 0 & datafull$z > 0),"fourth_dim"]
express_tran2 = datafull [which(datafull$y > 0),"fourth_dim"]
# Combine the vectors into a data frame

titles = c("No\nDuplications", "Ordinary\nDuplications", "AMY2B to\nAMY1", "AMY2B to\nAMY2A")

data <- data.frame(
  value = c(express_null, express_amy2b, express_tran1, express_tran2),
  category = factor(rep(titles, 
                        times = c(length(express_null), length(express_amy2b), length(express_tran1),length(express_tran2))))
)

data$category <- factor(data$category, levels = titles)


# Create the box plot with scatter points overlayed and showing the mean as a line
ggplot(data, aes(x = category, y = value)) +
  geom_boxplot(outlier.shape = NA) +  # Create the box plot without showing outliers
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +  # Add scatter points with some jitter to avoid overlap
  stat_summary(fun = mean, geom = "crossbar", width = 0.75, fatten = 2, color = "red", show.legend = FALSE) +  # Add mean line
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 15),  # Increase the size of the category text
    axis.title.y = element_text(size = 15)   # Optionally increase the size of the y-axis text
  ) + 
  labs(title = "", x = "", y = "Expression level")



















library(ggplot2)
library(tidyr) 


coeffs <- coef(summary(model))

# Extracting means (coefficients) and SE for our predictor variables
means <- coeffs[1:3, 1]
se <- coeffs[1:3, 2]
names <- rownames(coeffs[1:3, ])

# Create a dataframe for plotting
plot_data <- data.frame(names, means, se)


# Plot using ggplot2
library(ggplot2)
ggplot(plot_data, aes(x=names, y=means)) +
	geom_point(size=3) +  # Show the means as points
	geom_errorbar(aes(ymin=means-se, ymax=means+se), width=0.2) +  # Error bars for SE
	labs(title="", x="", y="Expression Level") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for clarity
	
	
library(dplyr)
df_long <- df %>%
  select(AMY_group1_0, AMY_group1_1, AMY_group1_2, express) %>%  # Selecting required columns
  pivot_longer(cols = -express, names_to = "Group", values_to = "Value")
  
  
  df_long <- df_long %>%
	mutate(express = log(express))
	






