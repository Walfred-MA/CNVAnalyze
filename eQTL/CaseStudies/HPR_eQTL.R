table <- read.csv("~/Documents/Marklab/eQTL/HPR_GETx_eQTL.tsv",sep = '\t', row.names=1)
table$Pop = NULL

table = data.frame(table)

expr_text = as.formula ( paste(c("HP ~ " , colnames(table)[1:(length(colnames(table))-2)] , " 0 "), collapse = " + " ) )
expr_text2 = as.formula ( paste(c("HPR ~ " , colnames(table)[1:(length(colnames(table))-2)] , " 0 "), collapse = " + " ) )
model1 = lm(expr_text , data = table)
model2 = lm(expr_text2 , data = table)

summary(model1)
summary(model2)

coef1 = coef(model1)
coef2 = coef(model2)

csums = colSums(table[,c(1:(ncol(table)-2))])

sumtable = cbind(coef1,coef2,csums )

new_rows <- data.frame(coef1 = rep(0, 3), coef2 = rep(0, 3), csums = rep(0, 3))

# Set row names for the new rows
rownames(new_rows) <- c("HP_group71_8", "HP_group71_9", "HP_group71_10")

sumtable  = rbind(sumtable [1:(9-1), ], new_rows, sumtable [9:nrow(sumtable), ])

cor.test(table$HP, table$HPR)
cor.test(coef1, coef2)

df = sumtable
df$Group = rownames(df)


usevalues1 = df$coef1[which(df$csums < 5)]
usevalues2 = df$coef2[which(df$csums < 5)]
df$coef1_norm <- (df$coef1 - min(usevalues1, na.rm = TRUE)) / (max(usevalues1, na.rm = TRUE) - min(usevalues1, na.rm = TRUE))
df$coef2_norm <- (df$coef2 - min(usevalues2, na.rm = TRUE)) / (max(usevalues2, na.rm = TRUE) - min(usevalues2, na.rm = TRUE))

# Set the normalized values to NA where 'csums' is less than 5
df$coef1_norm[df$csums < 5] <- NA
df$coef2_norm[df$csums < 5] <- NA

# Melt the normalized data frame to long format for ggplot
df$Group <- factor(df$Group, levels = df$Group)

# Melt the normalized data frame to long format for ggplot
df_long <- melt(df, id.vars = "Group", measure.vars = c("coef1_norm", "coef2_norm"))


# Display the melted data frame for verification

# Create the heatmap with a black frame around each cell and customized color legend
ggplot(df_long, aes(x = variable, y = Group, fill = value)) +
  geom_tile(color = "black") +  # Create the heatmap tiles with a black border around each cell
  scale_fill_gradientn(
    colors = c("blue", "yellow", "orange", "red", "darkred"),  # Distinctive multi-color gradient
    na.value = "white", 
    name = "expression"  # Set the color legend title to "expression"
  ) +
  theme_minimal() +  # Clean theme
  theme(
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank(),  # Remove y-axis title
    axis.text.x = element_blank(),   # Remove x-axis text (ticks)
    axis.text.y = element_blank(),   # Remove y-axis text (ticks)
    axis.ticks = element_blank(),    # Remove axis ticks
    plot.title = element_blank(),    # Remove plot title
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),   # Remove border around the plot
    legend.title = element_text(size = 14, face = "bold"),  # Increase font size of the legend title and make it bold
    legend.text = element_text(size = 12)  # Increase font size of the legend labels
  )