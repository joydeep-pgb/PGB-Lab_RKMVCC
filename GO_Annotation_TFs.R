
library(tidyverse)

SYM <- read_tsv("AT_SYM.txt")
GO <- read_tsv("Late_Nodule_GO.txt")
GO_exploded <- GO %>%
  separate_rows(GENE, sep = ",")

GO_exploded <- inner_join(GO_exploded, SYM, by = "GENE") %>% select(1:9, SYM)

# Assuming your exploded table is named 'GO_exploded'
GO_collapsed <- GO_exploded %>%
  group_by(GO_acc, term_type, Term, queryitem, querytotal, bgitem, bgtotal, pvalue, FDR) %>%
  summarize(GENE = paste(SYM, collapse = ","), .groups = "drop")


write_tsv(GO_collapsed, "21_dpi_GO.txt")

