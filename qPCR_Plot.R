library(tidyverse)

# Sample data
data <- data.frame(
  lncRNA = c(2.67, 1.37, 2.11, 2.19, 1.99, 3.14, 2.67, 4.28, 3.92),
  mRNA = c(0.88, 1.34, 1.72, 1.87, 3.04, 2.63, 5.32, 4.73, 4.18),
  group = factor(c("D", "D", "D", "S", "S", "S", "DS", "DS", "DS"), 
                 levels = unique(c("D", "S", "DS")))
)

# Group by and summarise
result <- data %>% group_by(group) %>% summarise(
  lncRNA.FC = mean(lncRNA), 
  mRNA.FC = mean(mRNA),
  lncRNA.SE = sd(lncRNA)/sqrt(length(lncRNA)), 
  mRNA.SE = sd(mRNA)/sqrt(length(mRNA))
)

print(result)

## For possition doge
library(reshape2)
data_melt <- melt(result, id.vars = "group", measure.vars = c("lncRNA.FC", "mRNA.FC"),
                  variable.name = "Type", value.name = "FoldChange")
# Add the standard error column
data_melt$SE <- c(result$lncRNA.SE, result$mRNA.SE)


# Reshape data into long format
data_long <- data.frame(
  group = rep(result$group, 2),
  FoldChange = c(result$lncRNA.FC, result$mRNA.FC),
  SE = c(result$lncRNA.SE, result$mRNA.SE),
  Type = rep(c('lncRNA', 'mRNA'), each = nrow(result))
)

plot <- ggplot(data_long, aes(x = group, y = FoldChange)) +
  geom_bar(stat = "identity", aes(fill = Type), position = position_dodge(width = 0.7), width = 0.5, colour = "black") +
  geom_errorbar(aes(ymin = FoldChange - SE, ymax = FoldChange + SE), width = 0.15, position = position_dodge(width = 0.7)) +
  facet_wrap(~ Type, scales = "fixed") +  # Use fixed scales to align both facets
  labs(x = "", y = "Relative Fold Change", title = "") +
  theme_classic() +
  scale_fill_manual(values = c('lncRNA' = 'skyblue', 'mRNA' = 'orange')) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2), expand = c(0, 0)) + 
  theme(axis.line = element_line(size = 0.8, colour = "gray60"),
        axis.ticks = element_line(size = 0.8, colour = "gray60"),
        axis.text.x = element_text(size = 15, colour = "black"), 
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"), 
        axis.ticks.length = unit(0.2, "cm"), 
        legend.position = "none")
# Save the plot
ggsave("facet_wrap_plot.png", plot = plot, dpi = 600, width = 6, height = 5)
