file2file = '~/Documents/MarkLab/Ctyper/Locityper/locityper-benchmarking/QV/full_illumina.csv_hap.csv' 

file1file = "~/Documents/MarkLab/Ctyper/Locityper/nonleaves_lociresults_273CMR_withname.txt"
# Load the data
file1 <- read.table(file1file, header =F, sep = " ")  # Replace with your file1 path
file2 <- read.table(file2file, header = T, sep = "\t")   # Replace with your file2 path



file1[,3] = as.numeric(file1[, 3])

# Calculate -log2 values, setting log2(0) to 100
# Filter out rows where file1[, 3] is -1
file1_filtered <- file1[file1[, 3] != -1, ]
file1_value = file1_filtered[,3]
# Calculate -log2 values for the filtered data, setting log2(0) to 100
log2_col1 <- -10*log10(ifelse(file1_value == 0, NA, file1_value))
log2_col1[is.na(log2_col1)] <- 100



file2_value = as.numeric(file2$div )
file2_value = file2_value[!is.na(file2_value)]
log2_col2 <- -10*log10(ifelse(file2_value== 0, NA, file2_value))
log2_col2[is.na(log2_col2)] <- 100


# Combine the data for plotting
combined_data <- data.frame(
  value = c(log2_col1, log2_col2),
  source = c(rep("File1 Column 3", length(log2_col1)),
             rep("File2 Div Column", length(log2_col2)))
)



# Create normalized histogram data for file1
hist1 <- hist(
  0.5^log2_col1, 
  breaks = 50,  # Number of bins
  plot = FALSE  # Do not plot yet
)
hist1$density <- hist1$counts / sum(hist1$counts)  # Normalize by total counts

# Create normalized histogram data for file2
hist2 <- hist(
  0.5^log2_col2, 
  breaks = 50,  # Same number of bins
  plot = FALSE  # Do not plot yet
)
hist2$density <- hist2$counts / sum(hist2$counts)  # Normalize by total counts

# Plot the first normalized histogram
plot(
  hist1$mids, hist1$density, 
  type = "h",  # Histogram-like line
  lwd = 2,  # Line width
  col = rgb(0, 0, 1, 0.7),  # Blue color
  main = "Normalized Combined Histograms of -log2 Values",
  xlab = "-log2 Values",
  ylab = "Density",
  xlim = range(c(hist1$breaks, hist2$breaks)),  # Ensure both histograms fit
  ylim = c(0, max(c(hist1$density, hist2$density))),  # Ensure enough space
  panel.first = grid()  # Add a grid to the background
)

plot(
  hist1$mids, hist1$density, 
  type = "h", 
  lwd = 2, 
  col = rgb(0, 0, 1, 0.7), 
  main = "Normalized Combined Histograms of -log2 Values",
  xlab = "-log2 Values",
  ylab = "Density",
  xlim = c(0.01, max(c(hist1$breaks, hist2$breaks))),  # Show only x > 0.01
  ylim = c(0, max(c(hist1$density, hist2$density))),
  panel.first = grid()
)





# Identify values equal to 0 in each dataset
zero_values_col1 <- log2_col1 >90
zero_values_col2 <- log2_col2 >90

# Keep original datasets for total count calculation
total_count1 <- length(log2_col1)
total_count2 <- length(log2_col2)

# Filter out values equal to 0 for plotting
plot_log2_col1 <- log2_col1[!zero_values_col1]
plot_log2_col2 <- log2_col2[!zero_values_col2]

# Create histogram data for filtered datasets
hist1 <- hist(
  plot_log2_col1, 
  breaks = 50, 
  plot = FALSE
)
hist2 <- hist(
  plot_log2_col2, 
  breaks = 50, 
  plot = FALSE
)

# Normalize counts to percentages based on total counts
hist1$percent <- (hist1$counts / total_count1) * 100
hist2$percent <- (hist2$counts / total_count2) * 100

# Plot the first histogram as percentages
plot(
  hist1$mids, hist1$percent, 
  type = "h", 
  lwd = 2, 
  col = rgb(0, 0, 1, 0.7), 
  main = "Percentage Histograms (Excluding Values = 0 from Plot)",
  xlab = "-log2 Values",
  ylab = "Percentage",
  xlim = range(c(hist1$breaks, hist2$breaks)), 
  ylim = c(0, max(c(hist1$percent, hist2$percent))),
  panel.first = grid()
)

# Overlay the second histogram as percentages
lines(
  hist2$mids, hist2$percent, 
  type = "h", 
  lwd = 2, 
  col = rgb(1, 0, 0, 0.7)
)

# Add a legend
legend(
  "topright",
  legend = c("File1 Column 3", "File2 Div Column"),
  col = c(rgb(0, 0, 1, 0.7), rgb(1, 0, 0, 0.7)),
  lwd = 2,
  bty = "n"
)

# Filter out values equal to 0 for plotting
plot_log2_col1 <- log2_col1[!zero_values_col1]
plot_log2_col2 <- log2_col2[!zero_values_col2]

# Create density estimates for smooth area shapes
density1 <- density(plot_log2_col1, adjust = 1)  # Adjust controls smoothness
density2 <- density(plot_log2_col2, adjust = 1)

# Plot the first density as an area
plot(
  density1,
  type = "n", 
  main = "",
  xlab = "QV",         # <-- Change here (x-axis label)
  ylab = "Density",    # <-- Change here (y-axis label)
  xlim = range(c(density1$x, density2$x)),
  ylim = range(c(density1$y, density2$y))
)

# Polygon for Ctyper (blue)
polygon(
  density1,
  col = rgb(0, 0, 1, 0.5),  # Semi-transparent blue
  border = NA
)

# Polygon for Locityper (red)
polygon(
  density2,
  col = rgb(1, 0, 0, 0.5),  # Semi-transparent red
  border = NA
)

# Add a legend
legend(
  "topright",
  legend = c("Ctyper", "Locityper"),
  fill = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)),
  border = NA,
  bty = "n"
)


length(which(file2$div==0))/length(file2$div)
length(which(file1[(file1[,3] >= 0),3]==0))/length(file1[(file1[,3] >= 0),3])

mean(file1_value)
mean(file2_value)
