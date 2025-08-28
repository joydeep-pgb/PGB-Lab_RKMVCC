library(tidyverse)
library(WGCNA)
## Interaction data
DF_UR <- read_csv("DF_UR.interaction.txt")
DF_RI <- read_csv("DF_RI.interaction.txt")
UR_RI <- read_csv("UR_RI.interaction.txt")

## Expression data
TPM <- read_tsv("TPM_adjested.norm.txt")

## -------------------------- DF_UR Correlation data ------------------------------
DF_UR <- DF_UR %>% filter(`Hybridization Energy` < -20)
DF.UR.lnc <- TPM %>% filter(transcript_id %in% DF_UR$`Query name`) %>% column_to_rownames("transcript_id")
DF.UR.tar <- TPM %>% filter(transcript_id %in% DF_UR$`Target name`) %>% column_to_rownames("transcript_id")

# Define numbers of genes and samples
nSamples <- ncol(DF.UR.lnc)
## Correlation
DF.UR.cor <- cor(t(DF.UR.lnc), t(DF.UR.tar), use = "p")
DF.UR.pval <- corPvalueStudent(DF.UR.cor, nSamples)
# Assuming EF.cor and EF.pval are your data frames
cor.data <- as.data.frame(DF.UR.cor)
pval.data <- as.data.frame(DF.UR.pval)
# Reshape the correlation matrix into long format
DF.UR.cor_long <- cor.data %>%
  rownames_to_column(var = "Module") %>%
  pivot_longer(cols = -Module, names_to = "Trait", values_to = "Correlation")

# Reshape the p-value matrix into long format
DF.UR.pval_long <- pval.data %>%
  rownames_to_column(var = "Module") %>%
  pivot_longer(cols = -Module, names_to = "Trait", values_to = "Pvalue")

# Merge the two long format data frames based on Module and Trait
DF_UR.cor <- inner_join(DF.UR.cor_long, DF.UR.pval_long, by = c("Module", "Trait"))

DF_UR.cor <- DF_UR.cor %>% filter(Correlation >= 0.55 | Correlation <= -0.55) %>% rename(lncRNA = Module, Target = Trait)

write_tsv(DF_UR.cor %>% filter(Correlation >= 0.55 | Correlation <= -0.55) %>% rename(lncRNA = Module, Target = Trait), "DF_UR.cor.txt")

## ----------------------------- DF_RI Correlation data ----------------------------------
DF_RI <- DF_RI %>% filter(`Hybridization Energy` < -20)

## DF_UR Correlation data
DF.RI.lnc <- TPM %>% filter(transcript_id %in% DF_RI$`Query name`) %>% column_to_rownames("transcript_id")
DF.RI.tar <- TPM %>% filter(transcript_id %in% DF_RI$`Target name`) %>% column_to_rownames("transcript_id")

# Define numbers of genes and samples
nSamples <- ncol(DF.RI.lnc)
## Correlation
DF.RI.cor <- cor(t(DF.RI.lnc), t(DF.RI.tar), use = "p")
DF.RI.pval <- corPvalueStudent(DF.RI.cor, nSamples)
# Assuming EF.cor and EF.pval are your data frames
cor.data <- as.data.frame(DF.RI.cor)
pval.data <- as.data.frame(DF.RI.pval)
# Reshape the correlation matrix into long format
DF.RI.cor_long <- cor.data %>%
  rownames_to_column(var = "Module") %>%
  pivot_longer(cols = -Module, names_to = "Trait", values_to = "Correlation")

# Reshape the p-value matrix into long format
DF.RI.pval_long <- pval.data %>%
  rownames_to_column(var = "Module") %>%
  pivot_longer(cols = -Module, names_to = "Trait", values_to = "Pvalue")

# Merge the two long format data frames based on Module and Trait
DF_RI.cor <- inner_join(DF.RI.cor_long, DF.RI.pval_long, by = c("Module", "Trait"))

DF_RI.cor <- DF_RI.cor %>% filter(Correlation >= 0.65 | Correlation <= -0.65) %>% rename(lncRNA = Module, Target = Trait)

write_tsv(DF_RI.cor %>% filter(Correlation >= 0.65 | Correlation <= -0.65) %>% rename(lncRNA = Module, Target = Trait), "DF_RI.cor.txt")

## ----------------------------- UR_RI Correlation data ----------------------------
UR_RI <- UR_RI %>% filter(`Hybridization Energy` < -20)

## DF_UR Correlation data
UR.RI.lnc <- TPM %>% filter(transcript_id %in% UR_RI$`Query name`) %>% column_to_rownames("transcript_id")
UR.RI.tar <- TPM %>% filter(transcript_id %in% UR_RI$`Target name`) %>% column_to_rownames("transcript_id")

# Define numbers of genes and samples
nSamples <- ncol(UR.RI.lnc)
## Correlation
UR.RI.cor <- cor(t(UR.RI.lnc), t(UR.RI.tar), use = "p")
UR.RI.pval <- corPvalueStudent(UR.RI.cor, nSamples)
# Assuming EF.cor and EF.pval are your data frames
cor.data <- as.data.frame(UR.RI.cor)
pval.data <- as.data.frame(UR.RI.pval)
# Reshape the correlation matrix into long format
UR.RI.cor_long <- cor.data %>%
  rownames_to_column(var = "Module") %>%
  pivot_longer(cols = -Module, names_to = "Trait", values_to = "Correlation")

# Reshape the p-value matrix into long format
UR.RI.pval_long <- pval.data %>%
  rownames_to_column(var = "Module") %>%
  pivot_longer(cols = -Module, names_to = "Trait", values_to = "Pvalue")

# Merge the two long format data frames based on Module and Trait
UR_RI.cor <- inner_join(UR.RI.cor_long, UR.RI.pval_long, by = c("Module", "Trait"))

UR_RI.cor <- UR_RI.cor %>% filter(Correlation >= 0.65 | Correlation <= -0.65) %>% rename(lncRNA = Module, Target = Trait)

write_tsv(UR_RI.cor %>% filter(Correlation >= 0.65 | Correlation <= -0.65) %>% rename(lncRNA = Module, Target = Trait), "UR_RI.cor.txt")
