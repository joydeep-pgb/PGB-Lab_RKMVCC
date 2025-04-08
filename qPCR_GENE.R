# Load necessary library
library(ggpubr)

# Data
data <- data.frame(
  lncRNA = c(1.4, 2.3, 1.5, 0.7, 0.9, 1.1, 1.6, 0.2, 1.1),
  group = factor(c("LF", "LF", "LF", "GF", "GF", "GF", "RF", "RF", "RF"),
                 levels = unique(c("LF", "GF", "RF")))
)

# Box plot with stats
ggboxplot(data, x = "group", y = "lncRNA",
          color = "black",
          fill = "group",
          ylab = "lncRNA expression",
          xlab = "Group",
          ggtheme = labs_pubr()) +
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("LF", "GF"), c("GF", "RF"), c("LF", "RF")),
                     label = "p.signif") +
  geom_jitter(
    position = position_jitter(width = 0.2), 
    size = 5, 
    color = "black",
    fill = "white",
    shape = 21
  ) +
  scale_fill_manual(values = c('LF' = '#a5de03', 'GF' = '#d6d760', 'RF' = '#01796f'))
