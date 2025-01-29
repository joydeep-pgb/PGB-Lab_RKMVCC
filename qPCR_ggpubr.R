# Load the necessary libraries
library(ggplot2)
library(ggpubr)
library(tidyr)

# Create the data frame
data <- data.frame(
  lncRNA = c(1.4, 2.3, 1.5, 0.7, 0.9, 1.1, 1.6, 0.2, 1.1),
  miRNA = c(0.7, 1.8, 1.4, 2.8, 2.7, 3.5, 4.2, 4.6, 5.1),
  mRNA = c(2.7, 3.9, 2.5, 0.8, 1.3, 0.5, 0.4, 0.2, 0.7),
  group = factor(c("LF", "LF", "LF", "GF", "GF", "GF", "RF", "RF", "RF"),
                 levels = unique(c("LF", "GF", "RF"))))

# Reshape the data to long format
data_long <- pivot_longer(data, cols = c(lncRNA, miRNA, mRNA), names_to = "type", values_to = "value")

ggbarplot(
  data_long, 
  x = "group", 
  y = "value", 
  fill = "type",
  add = "mean", 
  color = "black",
  width = 0.6,
  facet.by = "type", 
  scales = "fixed", # Keep consistent y-scales across facets (free_y/fixed)
  xlab = "", 
  ylab = "Relative Fold Change",
  title = "",
  size = 1 #Increase border size
) + 
  # Add error bars
  stat_summary(
    fun.data = mean_se, 
    geom = "errorbar", 
    width = 0.2,
    size = 1
  ) + 
  # Add jitter points
  geom_jitter(
    position = position_jitter(width = 0.2), 
    size = 5, 
    color = "black",
    fill = "white",
    shape = 21
  ) + 
  # Add statistical comparisons
  stat_compare_means(
    comparisons = list(c("LF", "GF"), c("GF", "RF"), c("LF", "RF")),
    label = "p.signif", 
    method = "t.test", 
    tip.length = 0.05,
    vjust = 0.5, 
    size = 7, 
    label.y = c(4.5, 5, 5.5)) + 
  # Apply custom theme and scales
  theme_classic() +
  scale_fill_manual(values = c('lncRNA' = '#a5de03', 'miRNA' = '#d6d760', 'mRNA' = '#01796f')) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2), expand = c(0, 0)) + 
  theme(
    axis.line = element_line(linewidth = 1, colour = "#4a4e69"),
    axis.ticks = element_line(linewidth = 1, colour = "#4a4e69"),
    axis.text.x = element_text(size = 25, colour = "black"), 
    axis.text.y = element_text(size = 25, colour = "black"),
    axis.title.x = element_text(size = 25, colour = "black"), 
    axis.title.y = element_text(size = 25, colour = "black"), 
    axis.ticks.length = unit(0.3, "cm"), 
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_blank())
