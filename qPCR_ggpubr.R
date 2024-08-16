# Load the necessary libraries
library(ggplot2)
library(ggpubr)
library(tidyr)

# Create the data frame
data <- data.frame(
  lncRNA = c(2.67, 1.37, 2.11, 2.19, 1.99, 3.14, 2.67, 4.28, 3.92),
  mRNA = c(0.88, 1.34, 1.72, 1.87, 3.04, 2.63, 5.32, 4.73, 4.18),
  group = factor(c("D", "D", "D", "S", "S", "S", "DS", "DS", "DS"), 
                 levels = unique(c("D", "S", "DS")))
)

# Reshape the data to long format
data_long <- pivot_longer(data, cols = c(lncRNA, mRNA), names_to = "type", values_to = "value")

# Create the faceted bar plot using ggpubr
ggbarplot(
  data_long, 
  x = "group", 
  y = "value", 
  fill = "type",
  add = "mean", 
  color = "black",
  width = 0.5,
  facet.by = "type", # Facet by the 'type' column
  scales = "fixed", # Allow different y-scales for each facet
  xlab = "", 
  ylab = "Relative Fold Change",
  title = ""
) + 
  # Add error bars
  stat_summary(
    fun.data = mean_se, 
    geom = "errorbar", 
    width = 0.2
  ) + 
  # Add statistical comparisons
  stat_compare_means(
    comparisons = list(c("D", "S"), c("D", "DS"), c("S", "DS")),
    label = "p.signif", # Show significance levels
    method = "t.test", 
    tip.length = 0.02,
    vjust = 0.5, 
    size = 5, 
    step.increase = 0.13) + 
  # Apply custom theme and scales
  theme_classic() +
  scale_fill_manual(values = c('lncRNA' = 'skyblue', 'mRNA' = 'orange')) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2), expand = c(0, 0)) + 
  theme(
    axis.line = element_line(linewidth = 0.8, colour = "gray60"),
    axis.ticks = element_line(linewidth = 0.8, colour = "gray60"),
    axis.text.x = element_text(size = 15, colour = "black"), 
    axis.text.y = element_text(size = 15, colour = "black"),
    axis.title.y = element_text(size = 15, colour = "black"), 
    axis.ticks.length = unit(0.2, "cm"), 
    legend.position = "none"
  )
