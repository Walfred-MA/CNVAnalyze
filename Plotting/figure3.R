
library(tidyr)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(grid)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(stringr)
library(scales)
 
t <- read.table("~/Documents/Marklab/benchfigure3/PangenomeAlleles_typefix.tsv_1KG_summary.txt_hwe.tsv", sep = '\t', header = F)



t$af <- (2*t$V5 + t$V6) / (2 * (t$V5 + t$V6 + t$V7))
t$het <- t$V6 / (t$V5 + t$V6 + t$V7)
pdf("~/Documents/Marklab/benchfigure3//HWE-plot.pdf")
#p <- ggplot(t, aes(af, het))
#p + geom_hex(bins=100)  +  scale_fill_gradient(low="gray", high="black", trans="log", breaks=c(1,10,100,1000,10000)) + labs(x="Allele frequency", y="Heterozygosity") + theme_bw(base_size=14)
#+ theme(legend.text=element_text(size=16)) + theme(axis.text=element_text(size=16))
p <- ggplot(t, aes(af, het)) +
  geom_hex(bins = 100) +
  scale_fill_gradient(low = "gray", high = "black", trans = "log", breaks = c(1, 10, 100, 1000, 10000)) +
  labs(x = "Allele frequency", y = "Heterozygosity") +
  theme_bw(base_size = 14) +
  theme(
    legend.text = element_text(size = 16),  # Increase legend text size
    legend.title = element_text(size = 16),  # Increase legend title size
    axis.text = element_text(size = 16),  # Increase axis text size
    axis.title = element_text(size = 16),  # Increase axis title size
    
  )
dev.off()
length(which(t$V3<0.05))



t <- read.table("~/Documents/Marklab/benchfigure3/PangenomeAlleles_typefix.tsv_1KG_summary.txt_hwe.tsv", sep = '\t', header = F)


# Calculate allele frequency and heterozygosity
t$af <- (2*t$V8 + t$V9) / (2 * (t$V8 + t$V9 + t$V10))
t$het <- t$V6 / (t$V8 + t$V9 + t$V10)
t$group <- ifelse(t$V3 < 0.05, "P<0.05", NA)

# Open a PDF device to save the plot
pdf("~/Documents/Marklab/benchfigure3/HWE-plot.pdf")

p <- ggplot(t, aes(x = af, y = het)) +
  geom_hex(data = subset(t, V3 >= 0.05), bins=100)  +  
  scale_fill_gradient(low="gray", high="black", trans="log", breaks=c(1,10,100,1000,10000)) + 
  #geom_point(aes(color = het), size = 1.5) + # Apply gradient based on heterozygosity or another metric
  geom_point(data = subset(t, !is.na(group)), aes(x = af, y = het, color = group), size = 0.5) + # Only plot points that meet the condition
  scale_color_manual(values = c("P<0.05" = "red"), name = "", labels = "P<0.05") + # Define colors and labels for the legend
  labs(x = "Allele frequency", y = "Heterozygosity", title = "Hardy weinberg Equilibrium of Pangenome-Alleles") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(angle = 0, hjust = 1, size = 20), # Explicitly increase font size here
        axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        plot.title = element_text(hjust = 0.5, size = 24)) + # Adjust y-axis labels size as needed
  coord_fixed(ratio = 1) +
  guides(color = guide_legend(title = "Significance"))
plot(p)
# Close the PDF device
dev.off()

pdf("~/Documents/Marklab/benchfigure3/CND-plot.pdf")
keys <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 32, 37)
values <- c(221333901, 21849790, 8628294, 238256, 137058, 35588, 22788, 14423, 8379, 6517, 5534, 4862, 4312, 3719, 2943, 2075, 1519, 778, 472, 194, 80, 60, 25, 13, 3, 2, 3, 1, 1)

# Combine into a data frame
data <- data.frame(keys, values)
data <- data[0:6,]
library(ggplot2)

ggplot(data, aes(x = factor(keys), y = values)) + 
  geom_bar(stat = "identity", fill = "blue") +
  labs(y = "Count", x = "Copy number", title = "Copy number distribution") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18), # Explicitly increase font size here
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 18)) + # Adjust y-axis labels size as needed
  scale_x_discrete(limits = as.character(0:6)) # Display only categories 0 to 10
dev.off()


t <- read.table("~/Documents/Marklab/benchfigure3/GTEX_summaryalelle_corr_1108.txt_hwe_nonpop.tsv", sep = '\t', header = F)
library(ggplot2)
library(ggExtra)
library(gridExtra)

# Calculate allele frequency and heterozygosity
t$af <- (2*t$V11 + t$V12) / (2 * (t$V11 + t$V12 + t$V13))
t$het <- t$V12 / (t$V11 + t$V12 + t$V13)

# Open a PDF device to save the plot
pdf("~/Documents/Marklab/benchfigure3/HWE-plot-gtex.pdf")

p <- ggplot(t, aes(x = af, y = het)) +
  geom_hex(data = subset(t, V3 >= 0.05), bins=100)  +  
  scale_fill_gradient(low="gray", high="black", trans="log", breaks=c(1,10,100,1000,10000)) + 
  #geom_point(aes(color = het), size = 1.5) + # Apply gradient based on heterozygosity or another metric
  geom_point(data = subset(t, V3 < 0.05), aes(x = af, y = het), color = "red", size = 0.5) + # Highlight points
  labs(x = "Allele frequency", y = "Heterozygosity") +
  coord_fixed(ratio = 1) +
  theme_bw(base_size = 14)
plot(p)
# Close the PDF device
dev.off()


pdf("~/Documents/Marklab/benchfigure3/Trio-plot.pdf")
clean_data <- read.table("~/Documents/Marklab/benchfigure3/PangenomeAlleles_typefix.tsv_1KG_summary.txt_trio.tsv", sep = '\t', header = F, na.strings = "")


clean_data$ColorFactor <- factor(ifelse(is.na(clean_data[,ncol(clean_data)]), NA, clean_data[,ncol(clean_data)]),
                                 levels = c(1, 2, 3),
                                 labels = c("chrX", "chrY", "Telomere"))

# Adding a column for actual colors to use, including grey for NAs
clean_data$Color <- ifelse(is.na(clean_data$ColorFactor), "grey",
                           ifelse(clean_data$ColorFactor == "X", "green",
                                  ifelse(clean_data$ColorFactor == "Y", "blue", "red")))
ggplot(clean_data, aes(x = 1:nrow(clean_data), y = clean_data[, 1])) +
  geom_point(aes(color = ColorFactor), na.rm = TRUE) +
  scale_color_manual(values = c("chrX" = "green", "chrY" = "blue", "Telomere" = "red"),
                     na.value = "grey") +
  labs(title = "Benchmarking based on trios", x = "Gene groups", y = "1 - F1-scores", color = "") + # Set color legend title to empty
  theme_minimal() +
  theme(legend.position = c(0.8, 0.8),
  		legend.text = element_text(size = 22),
        axis.text = element_text(angle = 0, hjust = 1, size = 22),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        plot.title = element_text(hjust = 0.5, size = 22))
        
dev.off()


lt <- read.csv("~/Documents/Marklab/benchfigure3/leaveoneout/nonleave_match.txt", header=F, sep = '\t')

lt <- lt[which(lt$V2 > 1000),]
nrow(lt)
length(which(lt$V1 <= lt$V3))
length(which(lt$V1 == 0))
lt$V5 =  pmin(1.0,lt$V1 / lt$V2)
mean(lt$V5)

lt <- read.csv("~/Documents/Marklab/benchfigure3/leaveoneout/leave_match.txt", header=F, sep = '\t')

lt <- lt[which(lt$V2 > 1000),]
nrow(lt)
length(which(lt$V1 <= lt$V3))
length(which(lt$V1 == 0))
# Filter data
lt$V5 =  pmin(1.0,lt$V1 / lt$V2)
mean(lt$V5)

library(ggplot2)
library(cowplot)

# Define the plots with increased axis label sizes
hist_top <- ggplot() + 
  geom_histogram(aes(lt$V2)) + 
  labs(x="", y="") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 18), 
        plot.margin = margin(10, 10, 10, 10))  # Adjust margins

hist_right <- ggplot() + 
  geom_histogram(aes(lt$V5)) + 
  labs(y="", x="") + 
  theme_bw() + 
  coord_flip() + 
  scale_y_continuous(breaks = seq(0, 500000, by = 100000)) +  # Double the y-axis scale
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 22), 
        axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1), 
        plot.margin = margin(10, 20, 10, 10))  # Adjust margins with more space on the right

nliP <- ggplot(lt, aes(V2, V5))
nliP3 <- nliP + 
  geom_hex(bins = 100) + 
  scale_fill_gradient(low = "gray", high = "black", trans = "log", breaks = c(1, 10, 100, 1000, 10000), guide = "none") + 
  labs(x = "", y = "") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 24), axis.text = element_text(size = 22), 
        axis.text.x = element_text(angle=90, vjust = 0.5), 
        plot.margin = margin(10, 10, 10, 10))  # Adjust margins

empty <- ggplot() + theme_void()

# Align the plots
aligned_top <- align_plots(hist_top, nliP3, align = "v", axis = "tblr")
aligned_right <- align_plots(hist_right, nliP3, align = "h", axis = "tblr")

# Arrange the plots with appropriate spacing
final_plot <- plot_grid(
  plot_grid(aligned_top[[1]], empty, ncol = 2, rel_widths = c(4, 2)),
  plot_grid(aligned_top[[2]], aligned_right[[1]], ncol = 2, rel_widths = c(4, 2)),
  ncol = 1, rel_heights = c(1.5, 4), align = "hv", axis = "tblr", label_size = 14
)

print(final_plot)







file = "~/Documents/Marklab/benchfigure3/HPRCto1KG_freq.tsv"
data  <- read.table(file , sep = ",")
data[] <- lapply(data, function(x) ifelse(is.nan(x), 1.0, x))
data <- data[!data[, 2] %in% c('EUR', 'SAS'), ]
data$ratio = data$V3/data$V4
data$log_pvalue = pmin(10,-log10(data$V10))


filterfile = "~/Documents/Marklab/benchfigure3/assemfiltersummary.txt"
filterdata  <- read.table(filterfile  , sep = " ", header = F, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
filterdata$V4 = (filterdata$V2-filterdata$V3)/filterdata$V3


filterdata = t(filterdata)
filterdata <- data.frame(filterdata,stringsAsFactors = FALSE)


names <- sapply(data$V1, function(x) {
  paste(strsplit(x, "_")[[1]][1:2], collapse = "_")
})
names= as.vector(gsub("-", ".", names))

filter_values = sum(as.vector(unlist(filterdata[unique(names)]["V2",])))/sum(as.vector(unlist(filterdata[unique(names)]["V3",])))
data$filter = unlist(filterdata[names]["V4",])



lm_model <- lm(V4 ~ V3, data = data)
summary_lm <- summary(lm_model)

# Extract coefficients and R-squared
intercept <- coef(lm_model)[1]
slope <- coef(lm_model)[2]
r_squared <- summary_lm$r.squared

# Create a data frame for the color palette
color_palette <- data.frame(
  V2 = c("AFR", "AMR", "EAS", "EUR", "SAS"),
  RGB = c("255,205,51", "255,61,61", "173,255,51", "100,235,255", "255,48,255")
)

# Convert RGB string to individual R, G, B components
color_palette$R <- as.integer(sapply(strsplit(as.character(color_palette$RGB), ","), `[`, 1)) / 255
color_palette$G <- as.integer(sapply(strsplit(as.character(color_palette$RGB), ","), `[`, 2)) / 255
color_palette$B <- as.integer(sapply(strsplit(as.character(color_palette$RGB), ","), `[`, 3)) / 255

# Match values in V2 to corresponding RGB values from the color palette
color_vector <- color_palette$RGB[match(data$V2, color_palette$V2)]

# Add color vector as a new column to the original data frame
data$RGB <- color_vector

ggplot(data, aes(x = V4, y = V3, color = factor(V2))) +
  geom_point(size = 2) +
  scale_color_manual(values = rgb(color_palette$R, color_palette$G, color_palette$B)) +
  labs(x = "", y = "", color = "Populations") +
  theme(
    legend.position = c(0.9, 0.3),  # Move legend slightly down
    legend.justification = c(1, 1),
    legend.text = element_text(size = 24),  # Increase legend text size
    legend.title = element_text(size = 24),  # Increase legend title size
    plot.background = element_rect(fill = "white", color = NA),  # White plot background
    panel.background = element_rect(fill = "white", color = NA),  # White panel background
    legend.background = element_rect(fill = "white", color = NA),  # White legend background
    panel.grid.major = element_line(color = "grey80"),  # Major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text = element_text(size = 24),  # Increase axis text size
    axis.title = element_text(size = 24)  # Increase axis title size
  ) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(breaks = seq(floor(min(data$V3)), ceiling(max(data$V3)), by = 1)) +
  scale_y_continuous(breaks = seq(floor(min(data$V4)), ceiling(max(data$V4)), by = 1)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Increase size of geom_smooth
  annotate("text", x =3, y = 8, label = paste("y =", round(intercept, 2), "+", round(slope, 2), "* x\n","R^2 =", round(r_squared, 2)), hjust = 0, vjust = 0, color = "black", size = 8)
  
lm_model <- lm(V4 ~ V3, data = data[which(data$filter == 0),])
summary_lm <- summary(lm_model)

# Extract coefficients and R-squared
intercept <- coef(lm_model)[1]
slope <- coef(lm_model)[2]
r_squared <- summary_lm$r.squared


length(which(data$V10<0.05 & data$filter==0))

length(which(data$filter==0))


file = "/Users/walfred/Documents/Marklab/benchfigure3/TtoTbench.csv"
data  <- read.table(file , sep = ",", row.names = 1)

ErrorRate = as.vector(as.numeric(data$V8[2:length(data$V8)]))

data$V5 = NULL
data$V8 = NULL
data <- as.data.frame(t(data))
# Convert from wide to long format
long_data <- gather(data, key = "Type", value = "Value", -Class)
# Set the factor levels for Type to the desired order
long_data$Type <- factor(long_data$Type, levels = c("=1", "=2", "<=5","<10", "<20",">=20"))
long_data$Value = as.numeric(long_data$Value)

type_order <- c("TP","Mistype","FP","FN","Misclassify")
long_data$Class <- factor(long_data$Class, levels = type_order)

totals <- long_data %>%
  filter(Type != "ErrorRate") %>%  # Exclude rows where Type is "ErrorRate"
  group_by(Type) %>%
  summarise(Total = sum(Value))

label_pos <- long_data %>%
  group_by(Type) %>%
  summarise(LabelPos = max(Value)) %>%
  ungroup()
  
# Choose a palette with distinguishable colors
palette <- brewer.pal(n = length(unique(long_data$Class))+4, name = "Set3")
palette <- palette[4:length(palette)]
# Plot a stacked bar plot with the new palette
ggplot(data = long_data, aes(x = Type, y = Value)) +
  geom_bar(aes(fill = Class), stat = "identity", position = "dodge", width = 0.4)  +
  geom_text(data = label_pos, aes(label = sprintf("%.2f%%", ErrorRate * 100), x = Type, y = LabelPos + 50), 
            color = "black", size = 8, vjust = 0) +  # Add labels
  scale_fill_manual(values = palette, labels = c("True Positive","Mistype","False Positive","False Negative", "Out of Reference")) +
  theme_minimal() +
  theme(axis.title = element_text(size = 24)) +
  labs(x = 'Number of Paralogs in CHM13', y = 'Counts', fill = '', size = 24) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(labels = c("=1", "=2", "3-5", "6-9", "10-19", ">=20")) +
  theme(axis.text = element_text(angle = 0, hjust = 0.5, size = 24),
        legend.text = element_text(size = 24), legend.title=element_text(size=24)) + 
  theme(legend.position = "right",legend.justification = c(1, 1),legend.margin = margin(t = 0, b = 0, l = -150, r = 10, unit = "pt"))






percent <- function(x, digits = 4) {
  paste0(round(x * 100, digits), "%")
}

leavesfile = "~/Documents/Marklab/benchfigure3/leaveoneout/leave_match.txt"
lt <- read.csv(leavesfile, header=F, sep = '\t')

lt <- lt[which(lt$V2 > 1000),]
leaves_v = unlist(lt$V1)/unlist(lt$V2)
closest_v = unlist(lt$V3)/unlist(lt$V2)

nonleavesfile = "~/Documents/Marklab/benchfigure3/leaveoneout/nonleave_match.txt"

lt2 <- read.csv(nonleavesfile, header=F, sep = '\t')
lt2 <- lt2[which(lt2$V2 > 1000),]
nonleaves_v = unlist(lt2$V1)/unlist(lt2$V2)

reffile = "~/Documents/Marklab/benchfigure3/leaveoneout/ref_match.txt"

lt3 <- read.csv(reffile, header=F, sep = '\t')
lt3 <- lt3[which(lt3$V2 > 1000 ),]
ref_v = unlist(lt3 $V1)/unlist(lt3$V2)



# Create a data frame with means and standard deviations
data <- data.frame(
  class = c("leaves-one-out genotype", "nonleaves-one-out genotype", "closest neightbor", "GRCh38"),
  Mean = c(mean(leaves_v), mean(nonleaves_v), mean(closest_v), mean(ref_v)),
  SD = c(sd(leaves_v), sd(nonleaves_v), sd(closest_v), sd(ref_v))
)


data$class <- factor(data$class, levels = c("leaves-one-out genotype","nonleaves-one-out genotype","closest neightbor", "GRCh38"))

label_position <- position_stack(vjust = 1)
new_labels <- c("leave-\none-out","nonleave-\none-out","closest\nneightbor", "GRCh38")
ggplot(data, aes(x = class, y = Mean, fill = class)) +
  geom_bar(stat = "identity") +
  labs(title = "", y = "Mean difference", x = "") +
  geom_text(aes(label = percent(Mean)), 
            color = "black", size = 8, vjust = 0) +  # Add labels
  scale_x_discrete(labels = new_labels) + 
  theme_minimal() +
  theme(axis.text = element_text(size = 24),  # Adjust x-axis label size
  			axis.text.y = element_text(size = 24),  # Adjust x-axis label size
  		axis.title = element_text(size = 20),
        legend.position = "none") +  # Clear legend
  guides(fill = FALSE)  # Remove legend guide  
 


# Assuming 'data' is your dataframe
timecost <- data.frame(
  coverage = c(20, 40, 60, 80, 100),
  time = c(4142.77, 5488.46, 6853.31, 7530.35,9208.72)
)

# Fit linear regression
lm_model <- lm(time ~ coverage, data = timecost)

# Extract coefficients of the linear model
intercept <- coef(lm_model)[1]
slope <- coef(lm_model)[2]

# Create scatter plot with linear trend line
p <- ggplot(timecost, aes(x = coverage, y = time)) +
  geom_point() +  # Add scatter points
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add linear trend line
  labs(title = "", x = "Coverage", y = "Time(s)") +
  theme(axis.text = element_text(size = 14),  # Adjust x-axis label size
  		axis.title = element_text(size = 20),
        legend.position = "none")   # Clear legend

# Calculate coordinates for placing the equation label
label_x <- max(timecost$coverage) * 0.9  # Adjust the multiplier as needed
label_y <- predict(lm_model, data.frame(coverage = label_x))

# Add equation label with formatted coefficients
p + geom_text(aes(x = label_x, y = label_y, 
                  label = paste("y =", round(slope, 2), "x +", round(intercept, 2)))) +
    geom_text(aes(label = round(time, 0)), vjust = -0.5, size = 5)








