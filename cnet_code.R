#load packages. If you install this tool locally, please Make sure that all of the following packages are installed.
library(shiny)
library(shinyjs)
library(shinyBS)
library(openxlsx)
library(gdata)
library(ggsci)
library(DT)
library(UpSetR)
library(glue)
library(ggplot2)
library(DOSE)
library(reshape2)
library(ggridges)
library(ComplexHeatmap)
library(clusterProfiler)
library(enrichplot)
library(circlize)
library(readxl)
library(tidyverse)


term <- read_excel("Term.xlsx", sheet = "Sheet2")
exp <- read_excel("Expression.xlsx")
fcvalue <- setNames(exp$Fold.Change, exp$IDs)

cnetplot.enrichResult(term, foldChange = fcvalue, 
                      layout = "kk", colorEdge = T, node_label = F, 
                      termcolor = "brown", palette = c("#4B0082", "#B3B3B3", "#CD5C5C"))

cnet(term, foldChange = fcvalue, 
                      layout = "kk", colorEdge = F, node_label = F, 
                      termcolor = "brown", palette = c("#ffd60a", "#339209"))
sessionInfo()
