install.packages("eulerr")

library(eulerr)

LFvsGF <- scan("LFvsGF.txt", what = "", sep = "\n")
GFvsRF <- scan("GFvsRipe.txt", what = "", sep = "\n")
LFvsRF <- scan("LFvsRipe.txt", what = "", sep = "\n")

venn_data <- list(
  LFvsGF = LFvsGF,
  GFvsRF = GFvsRF,
  LFvsRF = LFvsRF
)

venn_plot <- euler(venn_data)
plot(venn_plot, quantities = TRUE)
