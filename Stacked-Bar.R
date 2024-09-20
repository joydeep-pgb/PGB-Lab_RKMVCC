# Load necessary libraries
library(ggplot2)

# Create the dataset
data <- data.frame(
  Comparison = c("LFvsGF", "LFvsGF", "GFvsRF", "GFvsRF", "LFvsRF", "LFvsRF"),
  Direction = c("Up", "Down", "Up", "Down", "Up", "Down"),
  Count = c(139, 495, 426, 179, 371, 601)
)

# Set the order of the Comparison factor to match the order in the data
data$Comparison <- factor(data$Comparison, levels = c("LFvsGF", "GFvsRF", "LFvsRF"))
# Ensure "Down" appears below "Up" by ordering the levels of Direction
data$Direction <- factor(data$Direction, levels = c("Up", "Down"))

# Calculate the total for each Comparison
total_data <- data %>%
  group_by(Comparison) %>%
  summarise(Total = sum(Count))

# Create stacked bar plot with increased gap between comparisons
ggplot(data, aes(x = Comparison, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.3) +  # Reduce bar width
  scale_x_discrete(expand = expansion(mult = c(0.2, 0.01))) +  # Add gap between categories
  theme_classic() +
  labs(x = "", y = "Count", title = "") +
  scale_fill_manual(values = c("Up" = "#f0a967", "Down" = "#bbdd87")) +
  theme(axis.line = element_line(size = 1, colour = "#4a4e69"),
        axis.ticks = element_line(size = 1, colour = "#4a4e69"),
        axis.text.x = element_text(size = 25, colour = "black"), 
        axis.text.y = element_text(size = 25, colour = "black"),
        axis.title.x = element_text(size = 25, colour = "black"), 
        axis.title.y = element_text(size = 25, colour = "black"))

ggsave(width = 15, height = 6, dpi = 600,"UP-Down_Bar.png")
