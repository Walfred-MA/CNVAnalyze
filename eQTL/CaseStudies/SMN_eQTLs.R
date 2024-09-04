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

# Function to add size category
add_size_category <- function(data) {
  data %>%
    mutate(size_category = case_when(
      samplesize <= 2 ~ 1,
      samplesize <= 5 ~ 2,
      samplesize <= 10 ~ 4,
      samplesize <= 20 ~ 6,
      TRUE ~ 8
    ))
}

# Define a theme to center plot titles
center_title_theme <- theme(
  plot.title = element_text(hjust = 0.5)
)

library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot) # For get_legend

lines <- readLines("~/Documents/Marklab/eQTL/SMN1_SMN2_matrix.txt")
lines <- lines[lines != ""]
table = read.table(textConnection(lines), sep = ',', header =1, skip = 0, row.name = 1)
cols = colnames(table)
table = as.data.frame(t(table))
rownames(table) = cols
express = table[,length(colnames(table))]

#filter out deleted SMN2
table = table[which(table$SMN_group1_17_0==0),]
SMN1_ = "SMN_group1_0_0,SMN_group1_1_0,SMN_group1_2_0,SMN_group1_3_0,SMN_group1_9_0,SMN_group1_10_0,SMN_group1_14_0"
SMN2_ = "SMN_group1_4_0,SMN_group1_5_0,SMN_group1_6_0,SMN_group1_7_0,SMN_group1_8_0,SMN_group1_11_0,SMN_group1_12_0,SMN_group1_13_0"
SMN0_ = "SMN_group1_15_0,SMN_group1_16_0"

SMN1_ = strsplit(SMN1_, split = ",")[[1]]
SMN2_ = strsplit(SMN2_, split = ",")[[1]]
SMN0_ = strsplit(SMN0_, split = ",")[[1]]

table$express = table[,length(colnames(table))]
table$allSMN1 = rowSums2(table[, SMN1_])
table$allSMN2 = rowSums2(table[, SMN2_])
table$allSMN0 = rowSums2(table[, SMN0_])
table$allSMN12 = table$allSMN2 + table$allSMN1

model_formula <- as.formula("express ~ allSMN12  + allSMN0 + 0")

model <- lm(model_formula , data = table)

result = linearHypothesis(model, "allSMN12 - allSMN0" ,test='Chisq')
result


x = unlist(table$allSMN1)
y = unlist(table$allSMN2)
z = unlist(table$allSMN0)

fourth_dim = unlist(table$express)



# Combine the coordinates into a data frame
data <- data.frame(x, y, z, fourth_dim)

# Normalize the average fourth dimension values to a range [0, 1]
normalized_fourth_dim <- (fourth_dim - min(fourth_dim)) / (max(fourth_dim) - min(fourth_dim))
data <- data %>% distinct()
data$fourth_dim <- normalized_fourth_dim
# Create a color palette
colors <- colorRampPalette(c("blue", "red"))(256)

# Determine common x and y limits
common_lim <- range(c(data$x,data$y,data$z))


# Group by (x, y) and compute the average of the fourth dimension and sample size
avg_data_xy <- data %>%
  group_by(x, y) %>%
  summarise(normalized_fourth_dim = mean(fourth_dim), samplesize = n()) %>%
  ungroup() %>%
  add_size_category()

# Normalize the average fourth dimension values to a range [0, 1]

# Plot x-y plane without legend
plot_xy <- ggplot(avg_data_xy, aes(x = x, y = y, fill = normalized_fourth_dim, size = as.numeric(size_category))) +
  geom_point(shape = 21, alpha = 0.6) +  # Use shape 21 for filled circles
  scale_fill_gradientn(colors = colorRampPalette(c("blue", "red"))(256), guide = FALSE) +
  scale_size_continuous( guide = FALSE) +  # Adjust the range as needed and disable size guide
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_fixed() +  # Ensures x and y have the same scale
  theme_minimal() + center_title_theme + 
  labs(title = "Transcripts of SMN2 vs SMN1", x = "SMN1", y = "SMN2") +
  xlim(common_lim ) + ylim(common_lim )

# Group by (x, z) and compute the average of the fourth dimension and sample size
avg_data_xz <- data %>%
  group_by(x, z) %>%
  summarise(normalized_fourth_dim = mean(fourth_dim), samplesize = n()) %>%
  ungroup() %>%
  add_size_category()
# Normalize the average fourth dimension values to a range [0, 1]
# Plot x-z plane without legend
plot_xz <- ggplot(avg_data_xz, aes(x = x, y = z, fill = normalized_fourth_dim, size = as.numeric(size_category))) +
  geom_point(shape = 21, alpha = 0.6) +  # Use shape 21 for filled circles
  scale_fill_gradientn(colors = colorRampPalette(c("blue", "red"))(256), guide = FALSE) +
  scale_size_continuous( guide = FALSE) +  # Adjust the range as needed and disable size guide
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_fixed() +  # Ensures x and y have the same scale
  theme_minimal() + center_title_theme + 
  labs(title = "Transcripts of SMN-Converted vs SMN1", x = "SMN1", y = "SMN-Converted") +
  xlim(common_lim ) + ylim(common_lim )

# Group by (y, z) and compute the average of the fourth dimension and sample size
avg_data_yz <- data %>%
  group_by(y, z) %>%
  summarise(normalized_fourth_dim = mean(fourth_dim), samplesize = n()) %>%
  ungroup() %>%
  add_size_category()

# Plot y-z plane without legend
plot_yz <- ggplot(avg_data_yz, aes(x = y, y = z, fill = normalized_fourth_dim, size = as.numeric(size_category))) +
  geom_point(shape = 21, alpha = 0.6) +  # Use shape 21 for filled circles
  scale_fill_gradientn(colors = colorRampPalette(c("blue", "red"))(256), guide = FALSE) +
  scale_size_continuous( guide = FALSE) +  # Adjust the range as needed and disable size guide
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_fixed() +  # Ensures x and y have the same scale
  theme_minimal() + center_title_theme + 
  labs(title = "Transcripts of SMN-Converted vs SMN2", x = "SMN2", y = "SMN-Converted") +
  xlim(common_lim ) + ylim(common_lim )


# Extract the color legend from one of the plots
color_legend_plot <- ggplot(avg_data_xy, aes(x = x, y = y, fill = normalized_fourth_dim)) +
  geom_point(shape = 21, size = 3) +  # Use shape 21 for filled circles
  scale_fill_gradientn(colors = colorRampPalette(c("blue", "red"))(256)) +
  theme_minimal() +
  labs(fill = "Expression level")

# Extract the size legend from one of the plots
size_legend_plot <- ggplot(avg_data_xy, aes(x = x, y = y, size = factor(size_category))) +
  geom_point(shape = 21, fill = "grey", color = "black", alpha = 0.6) +  # Use shape 21 for filled circles
  scale_size_manual(values = c(`1` = 1,`2` = 2, `4` = 4, `6` = 6, `8` = 8, `10` = 10 ), 
  labels = c(`1` = "1-2",`2` = "3-5", `4` = '6-10', `6` = '11-20', `8` = '>20' )) + 
  #scale_size_manual(values = c( = 3, "6-10" = 4, "11-20" = 5, "21-50" = 6, ">50"= 7)) +
  
  theme_minimal() +
  labs(size = "Sample Size")
  
# Extract the legends
color_legend <- get_legend(color_legend_plot)
size_legend <- get_legend(size_legend_plot)









lines <- readLines("~/Documents/Marklab/eQTL/SMNSplice_matrix.txt")
lines <- lines[lines != ""]
table = read.table(textConnection(lines), sep = ',', header =1, skip = 0, row.name = 1)
cols = colnames(table)
table = as.data.frame(t(table))
rownames(table) = cols
express = table[,length(colnames(table))]

#filter out deleted SMN2
table = table[which(table$SMN_group1_17_0==0),]
SMN1_ = "SMN_group1_0_0,SMN_group1_1_0,SMN_group1_2_0,SMN_group1_3_0,SMN_group1_9_0,SMN_group1_10_0,SMN_group1_14_0"
SMN2_ = "SMN_group1_4_0,SMN_group1_5_0,SMN_group1_6_0,SMN_group1_7_0,SMN_group1_8_0,SMN_group1_11_0,SMN_group1_12_0,SMN_group1_13_0"
SMN0_ = "SMN_group1_15_0,SMN_group1_16_0"

SMN1_ = strsplit(SMN1_, split = ",")[[1]]
SMN2_ = strsplit(SMN2_, split = ",")[[1]]
SMN0_ = strsplit(SMN0_, split = ",")[[1]]

table$express = table[,length(colnames(table))]
table$allSMN1 = rowSums2(table[, SMN1_])
table$allSMN2 = rowSums2(table[, SMN2_])
table$allSMN0 = rowSums2(table[, SMN0_])
table$allSMN12 = table$allSMN2 + table$allSMN1

model_formula <- as.formula("express ~ allSMN1 + allSMN2  + allSMN0 + 0")

model <- lm(model_formula , data = table)

result = linearHypothesis(model, "allSMN2 - allSMN0" ,test='Chisq')
result


x = unlist(table$allSMN1)
y = unlist(table$allSMN2)
z = unlist(table$allSMN0)

fourth_dim = unlist(table$express)

# Combine the coordinates into a data frame
data <- data.frame(x, y, z, fourth_dim)

# Normalize the average fourth dimension values to a range [0, 1]
normalized_fourth_dim <- (fourth_dim - min(fourth_dim)) / (max(fourth_dim) - min(fourth_dim))
data <- data %>% distinct()
data$fourth_dim <- normalized_fourth_dim
# Create a color palette
colors <- colorRampPalette(c("blue", "red"))(256)

# Determine common x and y limits
common_lim <- range(c(data$x,data$y,data$z))

# Group by (y, z) and compute the average of the fourth dimension and sample size
avg_data_yz2 <- data %>%
  group_by(y, z) %>%
  summarise(normalized_fourth_dim = mean(fourth_dim), samplesize = n()) %>%
  ungroup() %>%
  add_size_category()

# Plot y-z plane without legend
plot_yz2 <- ggplot(avg_data_yz2, aes(x = y, y = z, fill = normalized_fourth_dim, size = as.numeric(size_category))) +
  geom_point(shape = 21, alpha = 0.6) +  # Use shape 21 for filled circles
  scale_fill_gradientn(colors = colorRampPalette(c("blue", "red"))(256), guide = FALSE) +
  scale_size_continuous( guide = FALSE) +  # Adjust the range as needed and disable size guide
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_fixed() +  # Ensures x and y have the same scale
  theme_minimal() + center_title_theme + 
  labs(title = "Exon7 Splicing of SMN-Converted vs SMN2", x = "SMN2", y = "SMN-Converted") +
  xlim(common_lim ) + ylim(common_lim )







# Arrange the four plots in a 2x2 grid
plots_grid <- arrangeGrob(
  plot_xy, plot_xz, 
  plot_yz, plot_yz2, 
  ncol = 2
)

# Combine the grid of plots with the legends
final_plot <- grid.arrange(
  arrangeGrob(plots_grid),
  arrangeGrob(color_legend, size_legend, ncol = 1),
  ncol = 2,
  widths = c(4, 1)
)













