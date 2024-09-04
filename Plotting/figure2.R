library(ggplot2)

thefile='~/Documents/Marklab/figure2/allPangenomeAlleles_locus.txt_num.txt'
data_raw=sort(unlist(read.table(thefile)))
data = log10(data_raw)
data  <- data.frame(index = 1:length(data), log10_sizes = data)

p <- ggplot(data, aes(x=index, y=log10_sizes)) +
    geom_point() +  # This adds the scatter plot points
    theme_minimal() +  # This applies a minimal theme to the plot
    labs(title="Sizes of pangenome-alleles", 
         x="Sorted Index", 
         y="Size") + 
    scale_y_continuous(breaks = log10(c(10, 100, 1000, 3000, 10000, 30000,100000,300000,1000000)), 
                       labels = c("10", "100", "1K", "3K","10K", "30K","100K","300K","1M")) +
    theme(text = element_text(size=12), 
    	  plot.title = element_text(hjust = 0.5),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
	theme(
    plot.title = element_text(hjust = 0.5, vjust = -2, size = 18),  # Center the title
    axis.title.x = element_text(size = 16, vjust = 2),  # Bold x-axis label
    axis.title.y = element_text(size = 16),  # Bold y-axis label
    axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5 , vjust = 5)
  )
  
  
  
ggsave("~/Documents/Marklab/figure2/figure2_allelesize.png", p, width = 5, height = 6, dpi = 300)



thefile='~/Documents/Marklab/figure2_typevsfreqs.csv'
data_raw=read.table(thefile,sep=',',header=F)

data = sort(unlist(data_raw[[1]]))
data  <- data.frame(index = 1:length(data), numbertypes = data)

p <- ggplot(data, aes(x=index, y= numbertypes)) +
    geom_point() +  # This adds the scatter plot points
    theme_minimal() +  # This applies a minimal theme to the plot
    labs(title="number of clade types in gene group", 
         x="Sorted Index", 
         y="number of clade types") + 
    theme(text = element_text(size=12), 
    	  plot.title = element_text(hjust = 0.5),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())

ggsave("~/Documents/Marklab/figure2_types.png", p, width = 10, height = 6, dpi = 300)




thefile='~/Documents/Marklab/figure2_typevsfreqs.csv'
data_raw=read.table(thefile,sep=',',header=F)

data = sort(unlist(data_raw[[2]]))
data  <- data.frame(index = 1:length(data), numberfreqs = data)

p <- ggplot(data, aes(x=index, y= numberfreqs)) +
    geom_point() +  # This adds the scatter plot points
    theme_minimal() +  # This applies a minimal theme to the plot
    labs(title="number of gene copies in gene group", 
         x="Sorted Index", 
         y="number of gene copies") + 
    theme(text = element_text(size=12), 
    	  plot.title = element_text(hjust = 0.5),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())

ggsave("~/Documents/Marklab/figure2_freqs.png", p, width = 10, height = 6, dpi = 300)



thefile='~/Documents/Marklab/figure2/allgenesizes.txt'
data=read.table(thefile,sep=' ',header=F)
data = log10(sort(unlist(data[[2]])))
   
# Create a histogram of the data
hist_obj <- hist(data, plot = FALSE)  # Generate histogram data but do not plot it yet

# Define the breaks for the x-axis
ticks <- hist_obj$breaks

# Determine which ticks are powers of ten
is_power_of_ten <- function(x) {
  x == as.integer(x) || x-as.integer(x) == 0.5
}
powers_of_ten <- sapply(ticks, is_power_of_ten)

# Calculate 10^x for those ticks only
labels <- ifelse(powers_of_ten, sapply(ticks, function(x) 10^x), NA)

# Format labels to remove scientific notation and assign NA where not a power of ten
formatted_labels <- ifelse(powers_of_ten, formatC(labels, format = "f", digits = 0, big.mark = ","), "")

thebreaks = log10(c(1000, 3000, 10000,30000, 100000,300000,1000000))
thelabels = c("1K", "3K", "10K", "30K", "100K", "300K", "1M")

pdf(file = "~/Documents/Marklab/figure2/allsizes.pdf") 

# Create the ggplot histogram
plot <- ggplot(data = data.frame(x = as.vector(data)), aes(x = x)) +
  geom_histogram(breaks = hist_obj$breaks, color = "black", fill = "skyblue") +
  labs(title = "Logistic distribution of pangenome allele sizes", x = "Allele Size", y = "Counts") +
  scale_x_continuous(breaks = thebreaks, labels = thelabels) +  # Setting custom x-axis labels
  theme_minimal() +  # Applying a minimal theme
  ylim(0, NA) +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = -2, size = 18),  # Center the title
    axis.title.x = element_text(size = 18, vjust = 2),  # Bold x-axis label
    axis.title.y = element_text(size = 18),  # Bold y-axis label
    axis.text.x = element_text(size = 18, angle = 0, hjust = 0.5 , vjust = 5),
    axis.text.y = element_text(size = 16, angle = 0, hjust = 0.5 , vjust = 5)
  )
  
# Print the plot
print(plot)
# Plot the histogram without the x-axis to add custom labels later
#hist(data, axes = FALSE, main = "Logistic distribution of pangenome allele sizes", xlab = "", col = "skyblue")
#hist(data, axes = FALSE, main = "Logistic distribution of pangenome allele sizes", xlab = formatted_labels, col = "skyblue")
#axis(side = 2)  # Add the y-axis
# Add customized x-axis with labels as integer values of 10^x where appropriate
#axis(side = 1, at = ticks, labels = formatted_labels, font = 2, las = 1)  # las=2 makes labels perpendicular to the axis

# Optionally, add a box around the plot
#box()
dev.off()


#sv composition

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(reshape2)
library(scales)



file = "/Users/walfred/Documents/Marklab/figure2/SVcalls.csv"
data  <- read.table(file , sep = ",", row.names = 1)
data <- data[, -which(colnames(data) == "V8")]
data <- as.data.frame(t(data))
data[,'Class'] = c("FullGene proximal Dup","Proximal Dup" ,"Exon deletion", "Deletion"  , "Insertion", "Small variants" )

# Convert from wide to long format
long_data <- gather(data, key = "Type", value = "Value", -Class)
# Set the factor levels for Type to the desired order
long_data$Type <- factor(long_data$Type, levels = c("Ref", "Alt", "Trans", "Novel"))
long_data$Value = as.numeric(long_data$Value)

totals <- long_data %>%
  group_by(Type) %>%
  summarise(Total = sum(Value))
  

# Choose a palette with distinguishable colors
palette <- brewer.pal(n = length(unique(long_data$Class))+4, name = "Set3")
palette <- palette[4:length(palette)]
# Plot a stacked bar plot with the new palette
ggplot(data = long_data, aes(x = Type, y = Value)) +
  geom_bar(aes(fill = Class), stat = "identity", position = "stack", width = 0.4) +
  geom_text(data = totals, aes(x = Type, label = scales::comma(Total), y = Total + 5000), 
            color = "black", size = 5, vjust = 0) +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 20)) +
  labs(x = NULL, y = 'Counts', fill = 'Variation Type', size = 20) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(labels = c("Reference\nClade", "Alternate\nClade", "Duplicated\nClade", "Novel\nClade")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),
        legend.text = element_text(size = 20), legend.title=element_text(size=20)) + 
  theme(legend.position = "right",legend.justification = c(1, 1),legend.margin = margin(t = 0, b = 0, l = -200, r = 10, unit = "pt"))
  
  
library(ggplot2)

file = "/Users/walfred/Documents/Marklab/figure2/genegroupinfo.txt"
data  <- read.table(file , sep = ",", col.names= c("group","Number of alleles","Number of allele-types", "Number of paralogs", "Number of CHM13 genes"))

# Categorizing 'Number of refgene' into bins
data$refgene_category <- cut(data$`Number.of.CHM13.genes`,
                             breaks = c(-Inf, 1, 2, 5, 10, 20, Inf),
                             labels = c("=1", "=2", "3-5", "6-10", "11-20", ">20"),
                             include.lowest = TRUE)

p <- ggplot(data, aes(x = `Number.of.allele.types`, y = `Number.of.paralogs`, size = `Number.of.alleles`, color = refgene_category)) +
  geom_point(alpha = 0.5) +  # Adjust transparency to see overlapping points
  theme_minimal() +
  labs(title = "All gene groups",
       x = "Number of allele-types",
       y = "Number of paralogs",
       color = "Number of CHM13 genes",
       size = "Number of alleles") +
  scale_size_continuous(range = c(0.1, 15), breaks = c(100, 200, 500, 1000, 2000, 4000, 6000)) +
  theme(plot.title = element_text(size=20, hjust = 0.5), axis.title = element_text(size = 18), 
  axis.text = element_text(size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 12)) +
  guides(color = guide_legend(title = "Number of CHM13 genes", override.aes = list(size = 6)),
         size = guide_legend(title = "Number of alleles"))
print(p)       

theme(plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),  # Increase legend title size
        legend.text = element_text(size = 12)) +  # Increase legend text size
  guides(
    color = guide_legend(title = "Number of CHM13 genes", override.aes = list(size = 6)),  # Increase symbol size in color legend
    size = guide_legend(title = "Number of alleles", override.aes = list(size = 6))  # Increase symbol size in size legend
  )
if (False)
{
	# Creating the scatter plot
p <- ggplot(data, aes(size = `Number.of.allele`, x = `Number.of.allele.types`, y = `Number.of.paralogs`, z =`Number of CHM13 genes` )) +
  geom_point(alpha = 0.5) +  # Adjust transparency to see overlapping points
  theme_minimal() +
  labs(title = "",
       color = "Number of alleles",
       x = "Number of allele-types",
       y = "Number of paralogs") +
  scale_size_continuous(range = c(0.1, 15), breaks = c(100,200,500,1000,2000,4000,6000)) +
  theme( plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 18))
# Display the plot

p <- p + guides(size = guide_legend(
  override.aes = list(alpha = 1),  # Make legend points fully opaque
  title = "Number of paralogs",
  ncol = 2,  # Single row for horizontal layout
  byrow = F,
  keywidth = 1.5,  # Adjust the width of keys
  keyheight = 1.5,  # Adjust the height of keys
  default.unit = "lines"
))
}










library(ggplot2)
library(gridExtra)
library(RColorBrewer)

# Assuming 'data' is your dataframe
data <- data.frame(
  Ref = c(561850, 3065, 655, 718, 652),
  Alt = c( 675977, 58360, 33163, 11766, 9432),
  Tran = c(29421, 1269, 1119, 269, 141),
  Novel = c(802, 3219, 8, 37, 49),
  Category = c("Small Variants", "Insertion", "Deletion", "Duplication", "Contraction")
)

# Prepare the color palette
palette <- brewer.pal(n = nrow(data), name = "Set3")

# Initialize list to store pie charts
pie_charts <- list()

# Iterate over each column except 'Category'
for (col in names(data)[-length(names(data))]) {
  # Prepare data for pie chart
  pie_data <- data.frame(Value = data[[col]], Category = data$Category)
  
  # Create pie chart
  p <- ggplot(pie_data, aes(x = "", y = Value, fill = Category)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = palette) +
    labs(title = col, x = NULL, y = NULL, fill = "Category") +
    theme_void() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  
  # Add pie chart to the list
  pie_charts[[col]] <- p
}

# Arrange the pie charts in a grid
do.call(grid.arrange, c(pie_charts, ncol = 2))



library(ggplot2)
library(tidyr)

#thefile = "~/Documents/Marklab/figure2/uniquekmers.csv"
thefile = "~/Documents/Marklab/figure2/group_summary.csv"
data_raw=data.frame(read.table(thefile,sep=',',header=F))
data_raw = data_raw[,c(1,2,3)]
colnames(data_raw) = c('name','all allele-types','allele-types with unique K-mers')
data_raw$name = NULL
data_raw$Index <- seq_len(nrow(data_raw))
# Convert data from wide to long format
data_long <- pivot_longer(data_raw, cols = -Index, names_to = "Variable", values_to = "Value")

# Create a mapping for custom color labels
custom_labels <- c(`all types` = "all within group", `allele-types with unique K-mers` = "allele-types with unique K-mers")

# Plot
ggplot(data_long, aes(x = Index, y = Value, color = Variable)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c(`all allele-types` = "blue", `allele-types with unique K-mers` = "red"),
                     labels = custom_labels) +
  labs(x = "Every Group", 
       y = "Counts", 
       color = "Number of allele-types", 
       title = "") + # Set the title of the plot
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0, size = 18),
        axis.text.y = element_text(angle = 0, hjust = 0, size = 18),
        axis.title = element_text(size = 18, hjust = 0.5),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.position = c(0.55, 0.8),
        plot.title = element_text(hjust = 0.5, size = 20))
        
        
        


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

thefile = "~/Documents/Marklab/figure2/all_dists.txt"
data_raw=data.frame(read.table(thefile,sep='\t',header=F))

# Assuming data_raw is your dataframe and has been loaded
# Rename columns V2 to V8 for clarity based on the provided names
colnames(data_raw)[2:8] <- c( 'reference', 'same\nallele-type', 'closest\nneighbor', '2nd\nclosest', '3rd\nclosest', '4th\nclosest', '5th\nclosest')

data_raw2 <- data_raw %>%
  mutate(across(colnames(data_raw)[2:8], ~na_if(., -1.0)))
  
 columns_order <- colnames(data_raw2)


# Reshape the filtered data to long format for plotting
long_data <- pivot_longer(data_raw2, cols = 2:8, names_to = "Metric", values_to = "Value")


# Transform 'Value' to its log10 equivalent, excluding non-positive values
long_data <- long_data %>%
  mutate(Value = ifelse(Value > 0, log10(Value), NA)) %>%
  drop_na() %>%
  mutate(Metric = factor(Metric, levels = columns_order))

# Plot all metrics in a single violin plot with log10-transformed values, preserving the order
ggplot(long_data, aes(x = Metric, y = Value, fill = Metric)) +
  geom_violin() +
  labs(title = "", x = "", y = "Logistic Differences") +  # Change y-axis label
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 18), 
  		axis.text.y = element_text(size = 18),  # Increase size of y-axis labels
        plot.title = element_text(hjust = 0.5, size = 15),
        axis.title = element_text(size = 18),
        legend.position = "none") +
  scale_y_continuous(labels = function(x) format(10^x, scientific = TRUE, trim = TRUE),  # Customize y-axis labels to show only one significant digit
                     breaks = seq(-5, 0, by = 1))  # Set breaks at every integer from -5 to 0





library(patchwork)


loadtable = function(file)
{
	table = read.csv(file, header = F)
	table[is.na(table)] = 0 
	thesum = colSums(table) - table[1,]
	table = rbind(table[1,], thesum)
	table =  t(table)
	table  = data.frame(table )
	colnames(table ) <- c("haplotype_num", "new_type_num")

	return (table )
}

 get_incre = function(thelist)
 {
 	incre = c()
 	incre[1] = thelist[1]
 	for(i in 2:(length(thelist))) {
  	incre[i] = thelist[i] - thelist[i-1]
	}
 	
 	return (incre)
 	
 }
 
NonAFR <- loadtable("~/Documents/Marklab/figure2/NonAFR_proj.csv") 

AFR <- loadtable("~/Documents/Marklab/figure2/AFR_proj.csv")  

# Initialize the first element to 0
NonAFR_2 <-loadtable("~/Documents/Marklab/figure2/NonAFR_proj_2.csv")  

NonAFR_incre = NonAFR_2 + AFR[1,2]
NonAFR_incre <- rbind(NonAFR_incre,AFR)
rownames(NonAFR_incre) <- rev(seq(length=nrow(NonAFR_incre)))
NonAFR_incre[,1]  =   rev(seq(length=nrow(NonAFR_incre)))


AFR_2  <-  loadtable("~/Documents/Marklab/figure2/AFR_proj_2.csv")  


AFR_incre = AFR_2 + NonAFR[1,2]
AFR_incre <- rbind(AFR_incre,NonAFR)
rownames(AFR_incre) <- rev(seq(length=nrow(AFR_incre)))
AFR_incre[,1] = rev(seq(length=nrow(AFR_incre)))
# Update the rest of the elements based on the difference from the previous element


NonAFR_incre$color_group <- ifelse(seq_len(nrow(NonAFR_incre)) > nrow(NonAFR), "African", "Non-African")
AFR_incre$color_group <- ifelse(seq_len(nrow(AFR_incre)) > nrow(AFR), "Non-African", "African")

AFR_sum = AFR
AFR_2_sum = AFR_2
AFR_2_sum$new_type_num = AFR_2_sum$new_type_num + NonAFR[1,2]
AFR_2_sum$haplotype_num <- rev(seq(length=nrow(AFR_2))+nrow(NonAFR))

NonAFR_sum = NonAFR
NonAFR_2_sum = NonAFR_2
NonAFR_2_sum$new_type_num = NonAFR_2_sum$new_type_num + AFR[1,2]
NonAFR_2_sum$haplotype_num <- rev(seq(length=nrow(NonAFR_2))+nrow(AFR))

p <- ggplot() + 
  geom_point(data=NonAFR_2_sum, aes(x=haplotype_num, y=new_type_num, color="Non-African"), size=1) + 
  geom_smooth(data=NonAFR_2_sum, aes(x=haplotype_num, y=new_type_num, color="Non-African"), method="lm", formula=y ~ poly(x, 21), size=1, se=FALSE) +
  geom_point(data=NonAFR_sum, aes(x=haplotype_num, y=new_type_num, color="Non-African"), size=0.5) + 
  geom_smooth(data=NonAFR_sum, aes(x=haplotype_num, y=new_type_num, color="Non-African"), method="lm", formula=y ~ poly(x, 21), size=0.5, se=FALSE) +
  # For African
  geom_point(data=AFR_2_sum, aes(x=haplotype_num, y=new_type_num, color="African"), size=1) + 
  geom_smooth(data=AFR_2_sum, aes(x=haplotype_num, y=new_type_num, color="African"), method="lm", formula=y ~ poly(x, 21), size=1, se=FALSE) +
  geom_point(data=AFR_sum, aes(x=haplotype_num, y=new_type_num, color="African"), size=0.5) + 
  geom_smooth(data=AFR_sum, aes(x=haplotype_num, y=new_type_num, color="African"), method="lm", formula=y ~ poly(x, 21), size=0.5, se=FALSE) +
  scale_color_manual(values = c("African" = "red", "Non-African" = "blue")) +
  # Add titles and labels
  labs(
    #title = "The number of total allele types as increase of sample size",
    x = "Number of haplotypes",
    y = "Number of allele clades",
    color = "Population"
  ) +
  # Center the plot title
  theme(plot.title = element_text(hjust = 0.5, size = 20)) + coord_cartesian(xlim = c(0, 250)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 18), 
  		axis.text.y = element_text(size = 18),  # Increase size of y-axis labels
        axis.title = element_text(size = 18), legend.title = element_text(size = 18), legend.text = element_text(size = 18),
        panel.background = element_rect(fill = "white", color = NA),  # White background
    panel.grid.major = element_line(color = "grey80"),  # Major grid lines
    panel.grid.minor = element_line(color = "grey90"),  # Minor grid lines
    plot.background = element_rect(fill = "white", color = NA),   # White plot background
    legend.title = element_text(size = 22),  # Increase size of legend title
    legend.text = element_text(size = 20),  # Increase size of legend text
    legend.position = c(0.8, 0.9)  # Move legend to the upper part of the plot
        )

# Display the plot



 
AFR_2$new_type_num = rev(get_incre(unlist(rev(AFR_2$new_type_num))))
NonAFR_2$new_type_num =  rev(get_incre(unlist(rev(NonAFR_2$new_type_num))))

AFR_2$haplotype_num = AFR_2$haplotype_num + nrow(NonAFR)
NonAFR_2$haplotype_num = NonAFR_2$haplotype_num + nrow(AFR)

AFR$new_type_num =  rev(get_incre(unlist(rev(AFR$new_type_num))))
NonAFR$new_type_num =  rev(get_incre(unlist(rev(NonAFR$new_type_num))))

q <- ggplot() + 
  # For Non-African
  geom_point(data=NonAFR_2, aes(x=haplotype_num, y=new_type_num, color="Non-African"), size=1) + 
  geom_smooth(data=NonAFR_2, aes(x=haplotype_num, y=new_type_num, color="Non-African"), method="lm", formula=y ~ poly(x, 21), size=1, se=FALSE) +
  #geom_point(data=NonAFR, aes(x=haplotype_num, y=new_type_num, color="Non-African"), size=0.5) + 
  #geom_smooth(data=NonAFR, aes(x=haplotype_num, y=new_type_num, color="Non-African"), method="lm", formula=y ~ poly(x, 21), size=0.5, se=FALSE) +
  # For African
  geom_point(data=AFR_2, aes(x=haplotype_num, y=new_type_num, color="African"), size=1) + 
  geom_smooth(data=AFR_2, aes(x=haplotype_num, y=new_type_num, color="African"), method="lm", formula=y ~ poly(x, 21), size=1, se=FALSE) +
  #geom_point(data=AFR, aes(x=haplotype_num, y=new_type_num, color="African"), size=0.5) + 
  #geom_smooth(data=AFR, aes(x=haplotype_num, y=new_type_num, color="African"), method="lm", formula=y ~ poly(x, 21), size=0.5, se=FALSE) +
  # Define custom colors for the legend
  scale_color_manual(values = c("Non-African" = "blue", "African" = "red")) +
  # Add titles and labels
  labs(
    #title = "The number of new discover allele types as increase of sample size",
    x = "Number of haplotypes",
    y = "Number of allele clades",
    color = "Population"
  ) +
  # Center the plot title
  theme(plot.title = element_text(hjust = 0.5, size = 20)) + coord_cartesian(xlim = c(0, 250))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 18), 
  		axis.text.y = element_text(size = 18),  # Increase size of y-axis labels
        axis.title = element_text(size = 18), legend.title = element_text(size = 18), legend.text = element_text(size = 18),
        panel.background = element_rect(fill = "white", color = NA),  # White background
    panel.grid.major = element_line(color = "grey80"),  # Major grid lines
    panel.grid.minor = element_line(color = "grey90"),  # Minor grid lines
    plot.background = element_rect(fill = "white", color = NA),   # White plot background
    legend.title = element_text(size = 22),  # Increase size of legend title
    legend.text = element_text(size = 20),  # Increase size of legend text
    legend.position = c(0.8, 0.9)  # Move legend to the upper part of the plot
        )

library(ggpubr)
# Display the plot
combined_plot = ggarrange(p, NULL,  q, ncol = 1, heights = c(1, 0.2, 1) )

# Display the combined plot
print(combined_plot)
















# Create a data frame with the values
data <- data.frame(
  Category = factor(c("TP", "FN", "FP"), levels = c("TP", "FN", "FP")),
  Count = c(3261, 215, 203)
)

# Create the bar plot
p <- ggplot(data, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 18),
    legend.position = "none"  # Hide legend as it's not needed
  ) +
  scale_fill_manual(values = c("TP" = "black", "FN" = "blue", "FP" = "red"))

# Print the plot
print(p)