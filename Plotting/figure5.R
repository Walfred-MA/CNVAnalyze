
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(grid)
library(gridExtra)
library(cowplot)
library(car)


qq_data <- function(p_values) {
  n <- length(p_values)
  expected <- -log10(ppoints(n))
  observed <- -log10(sort(p_values))
  data.frame(Expected = expected, Observed = observed)
}

file1="~/Documents/Marklab/eQTL/all_corr.txt"
file2="~/Documents/Marklab/eQTL/all_ftest.txt"
table1=read.table(file1, sep = ',', header = F)
table2=read.table(file2, sep = ',', header = F)
pvalues1 = as.numeric(unlist(table1$V3[which(table1$V3 != "NaN" & table1$V3 != 0.00)]))
pvalues2 = as.numeric(unlist(table2$V3[which(table2$V3 != "NaN" & table2$V3 != 0.00)]))

sign1 = table1[which(table1$V3 != "NaN" & table1$V3 != 0.00 & table1$V3 < 0.05/(length(pvalues1))),]
sign2 = table2[which(table2$V3 != "NaN" & table2$V3 != 0.00 & table2$V3 < 0.05/(length(pvalues2))),]

length(which(pvalues1 < 0.05/(length(pvalues1))))
length(which(pvalues2 < 0.05/(length(pvalues1))))
length(pvalues1)
length(pvalues2)

df1 <- qq_data(pvalues1)
df2 <- qq_data(pvalues2)
df1$Group <- "AggregateCN"
df2$Group <- "Allelespecific"
combined_df <- rbind(df1, df2)

ggplot(combined_df, aes(x = Expected, y = Observed, color = Group)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),  # Increase the size of the category text
    axis.title.y = element_text(size = 12), 
    legend.text = element_text(size = 12) ) + 
  labs(title = "", x = "Expected -log10(p)", y = "Observed -log10(p)") +
  scale_color_manual(name="",values = c("AggregateCN" = "blue", "Allelespecific" = "red"), labels = c("AggregateCN" = "Copy Number correlation\n(Spearman)\n", "Allelespecific" = "Allele specific vs AggregateCN\n(F-test)\n"))










qq_data <- function(p_values) {
  n <- length(p_values)
  expected <- -log10(ppoints(n))
  observed <- -log10(sort(p_values))
  data.frame(Expected = expected, Observed = observed)
}

file1 = "~/Documents/Marklab/eQTL/allpvalues_alt.txt"
file2 = "~/Documents/Marklab/eQTL/allpvalues_dup.txt"

data1 = read.table(file1,header = FALSE)
data2 = read.table(file2,header = FALSE)

pvalues1 = pmax(1e-300, as.numeric(unlist(data1[colnames(data1)[length(colnames(data1))]])))
pvalues2 = pmax(1e-300, as.numeric(unlist(data2[colnames(data2)[length(colnames(data2))]])))

df1 <- qq_data(pvalues1)
df2 <- qq_data(pvalues2)

# Add a column to identify the data sets
df1$Group <- "Alternative"
df2$Group <- "Duplicative"

# Combine the data frames
combined_df <- rbind(df1, df2)

ggplot(combined_df, aes(x = Expected, y = Observed, color = Group)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),  # Increase the size of the category text
    axis.title.y = element_text(size = 12), 
    legend.text = element_text(size = 12) ) + 
  labs(title = "", x = "Expected -log10(p)", y = "Observed -log10(p)") +
  scale_color_manual(name = "",values = c("Alternative" = "blue", "Duplicative" = "red"), labels = c("Alternative" = "Orthologs", "Duplicative" = "Paralogs")) +
  ylim(0,200)
    
length(which(pvalues2 < 0.05/(length(pvalues1)+length(pvalues2))))












file1 = "~/Documents/Marklab/eQTL/allpvalues_alt.txt"
file2 = "~/Documents/Marklab/eQTL/allpvalues_dup.txt"

data1 = read.table(file1,header = FALSE)
data2 = read.table(file2,header = FALSE)

pvalues1 = as.numeric(unlist(data1[colnames(data1)[length(colnames(data1))]]))
pvalues2 = as.numeric(unlist(data2[colnames(data2)[length(colnames(data2))]]))

data1 = data1[which(pvalues1 < 0.05/(length(pvalues1)+length(pvalues2))),]
data2 = data2[which(pvalues2 < 0.05/(length(pvalues1)+length(pvalues2))),]

pvalues1 = as.numeric(unlist(data1[colnames(data1)[length(colnames(data1))]]))
pvalues2 = as.numeric(unlist(data2[colnames(data2)[length(colnames(data2))]]))


ave1 = (data1$V8 * data1$V7 + data1$V5 * data1$V6)/(data1$V7 + data1$V6)
ave2 = (data2$V8 * data2$V7 + data2$V5 * data2$V6)/(data2$V7 + data2$V6)

afc1 = (data1$V8-data1$V6)/ave1 
afc2 = (data2$V8-data2$V6)/ave2

afc1[afc1 > 1] <- 1
afc1[afc1 < -1] <- -1
afc2[afc2 > 1] <- 1
afc2[afc2 < -1] <- -1

pvalues1[pvalues1 < 1e-100] <- -1
pvalues2[pvalues2 < 1e-100] <- -1


df <- data.frame(
  pvalue = c(pvalues1,pvalues2),
  effectsize = c(afc1, afc2),
  group = c(rep("Orthologs", length(afc1)), rep("Paralogs", length(afc2)))
)
df$log_pvalue <- -log10(df$pvalue)

dense_plot <- ggplot(df, aes(x = effectsize, color = group, y = ..density..)) +
  geom_freqpoly(bins = 51, size = 1) +
  theme_minimal() +
  labs(x = NULL, y = "Density") +
  scale_color_manual(name = "", values = c("blue", "red")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 15),
    legend.position = "none"
  )
  
# Prepare the data


# Create the volcano plot
volcano_plot = ggplot(df, aes(x = effectsize, y = log_pvalue, color = group)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  xlab("Effect Size") +
  ylab("-log10(p-value)") +
  scale_color_manual(name = NULL, values = c("blue", "red"))+
  ylim(0,101)+
  theme(
  axis.title.x = element_text(size = 15),
  axis.title.y = element_text(size = 15),
  legend.text = element_text(size = 15)
  )
  
  
# Combine the plots vertically
# Combine the plots vertically with adjusted height ratios
combined_plot <- dense_plot/ volcano_plot + plot_layout(heights = c(1, 3))


# Display the combined plot
print(combined_plot) 
  
# Perform the Wilcoxon rank-sum test
wilcox_test_result <- wilcox.test(abs(afc1), abs(afc2))

# Display the Wilcoxon test result
print(wilcox_test_result)









qq_data <- function(p_values) {
  n <- length(p_values)
  expected <- -log10(ppoints(n))
  observed <- -log10(sort(p_values))
  data.frame(Expected = expected, Observed = observed)
}


file1 = "~/Documents/Marklab/eQTL/expr_tissuesummary.txt_diff.txt_alt.txt"
file2 = "~/Documents/Marklab/eQTL/expr_tissuesummary.txt_diff.txt_dup.txt"

data1 = read.table(file1,header = FALSE)
data2 = read.table(file2,header = FALSE)

data1 = data1[which(data1$V11 != "NaN" & data1$V18 != "NaN"),]
data2 = data2[which(data2$V11 != "NaN" & data2$V18 != "NaN"),]

nrow_ori = 776902

pvalues1 = as.numeric(unlist(data1$V11)) + as.numeric(unlist(data1$V18)) -  as.numeric(unlist(data1$V11)) * as.numeric(unlist(data1$V18))
pvalues2 = as.numeric(unlist(data2$V11)) + as.numeric(unlist(data2$V18)) -  as.numeric(unlist(data2$V11)) * as.numeric(unlist(data2$V18))

pvalues1 = pmax(1e-300, pvalues1)
pvalues2 = pmax(1e-300, pvalues2)

df1 <- qq_data(pvalues1)
df2 <- qq_data(pvalues2)

# Add a column to identify the data sets
df1$Group <- "Alternative"
df2$Group <- "Duplicative"

# Combine the data frames
combined_df <- rbind(df1, df2)

ggplot(combined_df, aes(x = Expected, y = Observed, color = Group)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),  # Increase the size of the category text
    axis.title.y = element_text(size = 12), 
    legend.text = element_text(size = 12) ) + 
  labs(title = "", x = "Expected -log10(p)", y = "Observed -log10(p)") +
  scale_color_manual(name = "",values = c("Alternative" = "blue", "Duplicative" = "red"), labels = c("Alternative" = "Orthologs", "Duplicative" = "Paralogs")) +
  ylim(0,20)
  

sign1 = data1[(which(57*pvalues1 < 0.05/(length(pvalues1)+length(pvalues2)))),]
sign2 = data2[(which(57*pvalues2 < 0.05/(length(pvalues1)+length(pvalues2)))),]

sign1 = data1[which(pvalues1 < 0.05/(nrow_ori)),]
sign2 = data2[which(pvalues2 < 0.05/(nrow_ori)),]

sign1$group <- "Alternative"
sign2$group <- "Duplicative"

allsigns = rbind(sign1,sign2)

pvalues =  as.numeric(unlist(allsigns$V11)) + as.numeric(unlist(allsigns$V18)) - as.numeric(unlist(allsigns$V11)) * as.numeric(unlist(allsigns$V18))
effectsize = log(as.numeric(unlist(allsigns$V6))) - log(as.numeric(unlist(allsigns$V8))) + log(as.numeric(unlist(allsigns$V15))) - log(as.numeric(unlist(allsigns$V13)))

allsigns$log_pvalue <- pmin(20,-log10(pvalues))
allsigns$effectsize <- effectsize

# Create the volcano plot
ggplot(allsigns, aes(x = effectsize, y = log_pvalue, color = group)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  xlab("Effect Size") +
  ylab("-log10(p-value)") +
  scale_color_manual(name = NULL, values = c("blue", "red"))+
  ylim(0,21)+
  theme(
  axis.title.x = element_text(size = 15),
  axis.title.y = element_text(size = 15),
  legend.text = element_text(size = 15)
  )





file = "~/Documents/Marklab/eQTL/all_anova.txt_merge.txt"

data = read.table(file,header = FALSE,sep = ',')

data = data[which(data$V4 > 0.1 * data$V2 & data$V2 * data$V9  < 0.70*data$V4 & data$V3 < 46),]

var_total = data$V4
var_inst = (data$V2-1) * data$V9
var_random = data$V6
var_cnv = data$V7 
var_cnv_random = data$V8
var_explained = data$V5

snpnum = data$V11
panum = data$V3
mean_cnvs = data$V17

snp_total = pmax(0.00000001,data$V12)
snp_inst = (data$V10-1) * data$V16
snp_random = data$V14
snp_explained = data$V13

comb_total = pmax(0.00000001,data$V21)
comb_inst = (data$V10-1) * data$V25
comb_random = data$V23
comb_explained = data$V22

var_all_ratio = (var_explained - var_random)/(var_total - var_inst)
var_cnvs_ratio =  (var_cnv - var_cnv_random)/(var_total - var_inst)
var_type_ratio = (var_explained - var_cnv)/(var_total - var_inst)
snp_ratio = (snp_explained - snp_random)/(snp_total - snp_inst)
comb_ratio = ifelse(comb_explained > 0, (comb_explained - comb_random)/(var_total - var_inst),0)


var_all_ratio_uncorr = (var_explained)/(var_total - var_inst)
var_cnvs_ratio_uncorr =  (var_cnv )/(var_total - var_inst)
var_type_ratio_uncorr = (var_explained - var_cnv + (var_cnv_random+var_random)/2 )/(var_total - var_inst)
snp_ratio_uncorr = (snp_explained)/(snp_total - snp_inst)
comb_ratio_uncorr = ifelse(comb_explained > 0, (comb_explained)/(var_total - var_inst),0)

mean(comb_ratio[which(comb_ratio!=0)])
mean(comb_ratio_uncorr[which(comb_ratio_uncorr!=0)])
mean(snp_ratio[which(comb_ratio!=0)])

cor.test(snpnum,snp_ratio)
cor.test(snp_ratio,mean_cnvs)

table <- data.frame(
  x = 1:length(var_all_ratio),
  var_all_ratio = var_all_ratio,
  var_cnvs_ratio = var_cnvs_ratio,
  snp_ratio = snp_ratio,
  mean_cnv = mean_cnvs
)

# Categorize the mean_cnv values
table <- table %>%
  mutate(mean_cnv_category = case_when(
    mean_cnv < 0.01 ~ "< 0.01",
    mean_cnv >= 0.01 & mean_cnv < 0.1 ~ "0.01 - 0.1",
    mean_cnv >= 0.1 & mean_cnv < 0.5 ~ "0.1 - 0.5",
    mean_cnv >= 0.5 ~ ">0.5"
  ))

# Set the levels for the mean_cnv_category factor
table$mean_cnv_category <- factor(table$mean_cnv_category, levels = c("< 0.01", "0.01 - 0.1", "0.1 - 0.5", ">0.5"))

# Calculate the mean of each ratio for each category
mean_values <- table %>%
  group_by(mean_cnv_category) %>%
  summarise(
    mean_var_all_ratio = mean(var_all_ratio, na.rm = TRUE),
    mean_var_cnvs_ratio = mean(var_cnvs_ratio, na.rm = TRUE),
    mean_snp_ratio = mean(snp_ratio, na.rm = TRUE),
    mean_nozero_ratio = mean(snp_ratio[which(snp_ratio != 0)], na.rm = TRUE)
  )
  
# Reshape the data for plotting
mean_values_long <- mean_values %>%
  pivot_longer(cols = c(mean_var_all_ratio, mean_var_cnvs_ratio,  mean_snp_ratio, mean_nozero_ratio ), 
               names_to = "variable", 
               values_to = "value")

table_uncorr <- data.frame(
  x = 1:length(var_all_ratio_uncorr),
  var_all_ratio = var_all_ratio_uncorr,
  var_cnvs_ratio = var_cnvs_ratio_uncorr,
  snp_ratio = snp_ratio_uncorr,
  mean_cnv = mean_cnvs
)

# Categorize the mean_cnv values
table_uncorr <- table_uncorr %>%
  mutate(mean_cnv_category = case_when(
    mean_cnv < 0.01 ~ "< 0.01",
    mean_cnv >= 0.01 & mean_cnv < 0.1 ~ "0.01 - 0.1",
    mean_cnv >= 0.1 & mean_cnv < 0.5 ~ "0.1 - 0.5",
    mean_cnv >= 0.5 ~ ">0.5"
  ))

# Set the levels for the mean_cnv_category factor
table_uncorr$mean_cnv_category <- factor(table_uncorr$mean_cnv_category, levels = c("< 0.01", "0.01 - 0.1", "0.1 - 0.5", ">0.5"))

# Calculate the mean of each ratio for each category
mean_values_uncorr <- table_uncorr %>%
  group_by(mean_cnv_category) %>%
  summarise(
    mean_var_all_ratio = mean(var_all_ratio, na.rm = TRUE),
    mean_var_cnvs_ratio = mean(var_cnvs_ratio, na.rm = TRUE),
    #mean_var_type_ratio = mean(var_type_ratio, na.rm = TRUE),
    mean_snp_ratio = mean(snp_ratio, na.rm = TRUE),
    mean_nozero_ratio = mean(snp_ratio[which(snp_ratio != 0)], na.rm = TRUE)
  )
  
# Reshape the data for plotting
mean_values_uncorr_long <- mean_values_uncorr %>%
  pivot_longer(cols = c(mean_var_all_ratio, mean_var_cnvs_ratio,mean_snp_ratio, mean_nozero_ratio ), 
               names_to = "variable", 
               values_to = "uncorrected")                              

               
mean_values_both = mean_values_long
mean_values_both$uncorrected = mean_values_uncorr_long$uncorrected

# Update the variable names for better readability and set the order
mean_values_both$variable <- factor(mean_values_both$variable, 
                                    levels = c("mean_var_all_ratio", "mean_var_cnvs_ratio", "mean_snp_ratio", "mean_nozero_ratio"),
                                    labels = c("Allele specific copy numbers", "Aggregate copy numbers", "SNVs sites (All)", "SNVs sites (only eGenes)"))


 # Custom legend entry
custom_legend_data <- data.frame(
  mean_cnv_category = c(0,0),
  value = c(0,0),
  uncorrected = c(0,0),
  variable = factor(c("Corrected values", "Baseline by random genotypes"), levels = c(levels(mean_values_both$variable), "Corrected values", "Baseline by random genotypes")))

combined_data <- rbind(mean_values_both, custom_legend_data)


# Generate dynamic color mapping
color_mapping <- c(
  "Allele specific copy numbers" = "blue",
  "Aggregate copy numbers" = "green",
  "SNVs sites (All)" = "red",
  "SNVs sites (only eGenes)" = "orange"
)


# Create the bar plot with uncorrected values and a dotted area for the difference
p = ggplot(mean_values_both, aes(x = mean_cnv_category, fill = variable)) +
  geom_bar(aes(y = value), stat = "identity", position = "dodge") +
  geom_bar(aes(y = uncorrected), stat = "identity", position = "dodge", alpha = 0.5) +
  labs(
    title = "",
    x = "MAD of CNV",
    y = "Mean explained expression variances"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20), # Center title and make it larger
    axis.title = element_text(size = 18), # Increase axis titles size
    axis.text = element_text(size = 15), # Increase axis text size
    legend.title = element_text(size = 18), # Increase legend title size
    legend.text = element_text(size = 15) # Increase legend text size
  ) +
  scale_fill_manual(name="",values = color_mapping) +
  guides(fill = guide_legend(override.aes = list(alpha = c(rep(1, 4))))) 
 
 
# Create a dummy plot to extract the alpha legend
alpha_legend_plot <- ggplot() +
  geom_point(aes(x = 1, y = 1, fill = "Baseline subtracted"), shape = 15, size = 8, alpha = 1) +
  geom_point(aes(x = 1, y = 2, fill = "Baseline(Permutation)"), shape = 15, size = 8, alpha = 0.5) +
  scale_fill_manual(
    name = "",
    values = c("Baseline subtracted" = "blue", "Baseline(Permutation)" = "blue")
  ) +
  theme_void() +
  guides(fill = guide_legend(override.aes = list(alpha = c(1, 0.1)))) +
  theme(
    legend.position = c(0.48, 0.8),  # Adjust this to move the legend vertically
    legend.key.size = unit(1.5, "lines"),  # Increase the size of legend keys
    legend.text = element_text(size = 14)  # Increase the size of the legend text
  )
  
# Extract the legends
original_legend <- get_legend(p)
alpha_legend <- get_legend(alpha_legend_plot)

# Combine the original plot with the legends
combined_plot <- plot_grid(
  plot_grid(p + theme(legend.position = "none"), ncol = 1),
  plot_grid(original_legend, alpha_legend, ncol = 1, rel_heights = c(2, 1)),
  rel_widths = c(3, 1)
)

# Print the combined plot
print(combined_plot)


table <- data.frame(
  x = 1:length(var_all_ratio),
  var_all_ratio = sort ( pmin(1.0,  var_all_ratio),decreasing = T),
  var_cnvs_ratio = sort ( pmin(1.0,  var_cnvs_ratio),decreasing = T),
  snp_ratio = sort( pmin(1.0, snp_ratio),decreasing = T)
)

# Transform data into long format
data_long <- pivot_longer(table, cols = -x, names_to = "variable", values_to = "value")

# Ensure the variable factor has the correct order
data_long$variable <- factor(data_long$variable, levels = c("var_all_ratio", "var_cnvs_ratio",  "snp_ratio"))

ggplot(data_long, aes(x = x, y = value, color = variable)) +
  geom_line() +
  labs(
    title = "",
    x = "",
    y = "Explained expression variance"
  ) +
  theme_minimal(base_size = 15) + # Increase base text size
  theme(
    plot.title = element_text(hjust = 0.5, size = 20), # Center title and make it larger
    axis.title = element_text(size = 18), # Increase axis titles size
    axis.text = element_text(size = 15), # Increase axis text size
    legend.title = element_text(size = 18), # Increase legend title size
    legend.text = element_text(size = 15), # Increase legend text size
    legend.position = c(0.5, 0.5), # Adjust the legend position
    legend.justification = "left" # Justify the legend to the left
  ) +
  scale_color_manual(name = "",
    values = c(
      "var_all_ratio" = "blue",
      "var_cnvs_ratio" = "red",
      "snp_ratio" = "black"
    ),
    labels = c(
      "var_all_ratio" = "Allele specific copy numbers",
      "var_cnvs_ratio" = "Aggregate copy numbers",
      "snp_ratio" = "SNVs sites"
    )
  )












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






