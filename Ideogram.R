library(RIdeogram)

## Data loading
mango_karyotype <- read.table("karyotype.txt", sep = "\t", header = T, stringsAsFactors = F)
lncRNA_density <- read.table("lncRNA-gene_density.txt", sep = "\t", header = T, stringsAsFactors = F)
mRNA_density <- read.table("mRNA-gene_density.txt", sep = "\t", header = T, stringsAsFactors = F)
mRNA_3k <- read.table("mRNA_DEN.txt", sep = "\t", header = T, stringsAsFactors = F)

## Density extraction, can also be done directly from GTF file
gene_density <- GFFex(input = "DEL.Mango.GFF3.txt", karyotype = "karyotype.txt", feature = "gene", window = 1000000)


## Generate ideogram
ideogram(karyotype = mango_karyotype, overlaid = gene_density)

ideogram(karyotype = mango_karyotype, overlaid = lncRNA_density, label = mRNA_density, label_type = "heatmap", colorset1 = c("#a2d2ff", "#2c7fb8"), colorset2 = c("#efcbd0", "#ef233c"),  width = 170, Lx = 190, Ly = -2)

ideogram(karyotype = mango_karyotype, overlaid = lncRNA_density, label = mRNA_density, label_type = "heatmap", colorset1 = c("#aaf3ed", "#00afb9"), colorset2 = c("#ffdbb1", "#ff8800"),  width = 170, Lx = 190, Ly = -2)

svg2pdf("chromosome.svg", file = "Combined_3.pdf", dpi = 300)

