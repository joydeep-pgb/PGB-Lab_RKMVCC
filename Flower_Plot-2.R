## Custom ploting function
flower_plot <- function(sample, value, start, a, b,  
                        circle_col = rgb(255, 255, 255, max = 255),
                        circle_text_cex = 1.1, labels = labels) {
  par(bty = "n", ann = F, xaxt = "n", yaxt = "n", mar = c(1, 1, 1, 1))
  plot(c(0, 10), c(0, 10), type = "n", asp = 1)  # Set asp = 1 for equal aspect ratio
  n   <- length(sample)
  deg <- 360 / n
  
  # Adjust the ellipse colors using a gradient palette with transparency
  base_colors <- colorRampPalette(c('yellow', 'green', 'blue', 'magenta', 'red'))(n)
  gradient_colors <- sapply(base_colors, function(col) adjustcolor(col, alpha.f = 0.15))  # 0.5 transparency
  
  res <- lapply(1:n, function(t) {
    plotrix::draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180),
                          y = 5 + sin((start + deg * (t - 1)) * pi / 180),
                          col = gradient_colors[t],
                          border = gradient_colors[t],
                          a = a, b = b, angle = deg * (t - 1))
    text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
         value[t])
    
    if (deg * (t - 1) < 180 && deg * (t - 1) > 0) {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) - start,
           adj = 1,
           cex = circle_text_cex)
    } else {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) + start,
           adj = 0,
           cex = circle_text_cex)
    }
  })
  
  plotrix::draw.circle(x = 5, y = 5, r = 1.0, col = circle_col, border = circle_col)
  
  # Tune location by x and y.
  text(x = 4.7, y = 5, labels = labels)
}

# Call the flower_plot function with your data
flower_plot(
  c(
    "V. vinifera", "M. domestica", "A. chinensis", "J. regia", "C. clementina", "P. trichocarpa", "C. quinoa", "V. radiata", "P. persica", "M. truncatula", "C. avellana", "S. tuberosum", "I. triloba", "M. esculenta", "C. sativa", "C. capsularis", "G. raimondii", "B. rapa", "G. max", 
    "A. thaliana", "T. aestivum", "S. bicolor", "C. annuum", "O. sativa", 
    "Z. mays"),
  c(155, 150, 148, 141, 138, 137, 123, 122, 114, 112, 111, 105, 101, 99, 93, 87, 81, 80, 71, 70, 61, 39, 38, 22, 21), 
  90, 0.9, 2.0, 
  labels = "6315",
  circle_col = "#FFFFFF"
)

## More Perfect dimensions
flower_plot(
  c(
    "V. vinifera", "M. domestica", "A. chinensis", "J. regia", "C. clementina", "P. trichocarpa", "C. quinoa", "V. radiata", "P. persica", "M. truncatula", "C. avellana", "S. tuberosum", "I. triloba", "M. esculenta", "C. sativa", "C. capsularis", "G. raimondii", "B. rapa", "G. max", 
    "A. thaliana", "T. aestivum", "S. bicolor", "C. annuum", "O. sativa", 
    "Z. mays"),
  c(155, 150, 148, 141, 138, 137, 123, 122, 114, 112, 111, 105, 101, 99, 93, 87, 81, 80, 71, 70, 61, 39, 38, 22, 21), 
  90, 0.65, 1.75,  # Updated ellipse dimensions
  labels = "6315",
  circle_col = "#FFFFFF"
)
