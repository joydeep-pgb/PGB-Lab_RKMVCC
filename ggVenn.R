set.seed(20190708)
genes <- paste("gene",1:1000,sep="")
x <- list(
  A = sample(genes,300), 
  B = sample(genes,525), 
  C = sample(genes,440),
  D = sample(genes,350)
)

library(VennDiagram)

venn.diagram(x, filename = "venn-4-dimensions.png")

# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

# Three dimension Venn plot
display_venn(x[1:3])

# Further customization
display_venn(
  x[1:3],
  category.names = c("Set 1" , "Set 2 " , "Set 3"),
  # Circles
  lwd = 3,
  col = "white",
  fill = c("#4daf4a", "#fdc500", "#e6550d"), alpha = 0.6)
?venn.diagram()
