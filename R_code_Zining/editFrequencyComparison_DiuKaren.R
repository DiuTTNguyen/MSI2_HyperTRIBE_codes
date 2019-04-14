# Histogram of log2(fold change) of edited frequency of LSC vs LSK in shared targets
# Some modifications made by Diu and Karen
# Diu: Make sure you install all packages below with ie. install.packages("tidyverse")

library("ggplot2")
library("openxlsx")
library("stringr")
library("tidyverse")
library("VennDiagram")


# filenames - edits
#fname.edit.lsc <- file.path("hypertribe_data", "mouse_lsc_snp_counts_dedupped_significant.csv")
#fname.edit.lsk <- file.path("hypertribe_data", "mouse_lsk_snp_counts_dedupped_significant.csv")

# filenames - gene expression
#fname.gene.lsc <- file.path("hypertribe_data", "mouse_lsc_gene_expression.xlsx")
#fname.gene.lsk <- file.path("hypertribe_data", "mouse_lsk_gene_expression.xlsx")

# set working directory - Karen Chu
setwd( "/Users/chuk/Documents/MSI2-hyperTRIBE/data_Yuhueng")

# read CSV files
raw.edit.lsc <- read.csv("mouse_lsc_snp_counts_dedupped_significant.csv")
raw.edit.lsk <- read.csv("mouse_lsk_snp_counts_dedupped_significant.csv")

# read gene expression files
raw.gene.lsc <- openxlsx::read.xlsx("mouse_lsc_gene_expression.xlsx", sheet = 1)
raw.gene.lsk <- openxlsx::read.xlsx("mouse_lsk_gene_expression.xlsx", sheet = 1)

# figure out how many editing sites
# 1. create column that gives ID of each site (chromosome + position)
raw.edit.lsc$site.id <- paste(raw.edit.lsc$chr,  raw.edit.lsc$pos)
raw.edit.lsk$site.id <- paste(raw.edit.lsk$chr,  raw.edit.lsk$pos)

# define significance criteria for edits
p.thres <- 0.01 # cutoff p-value for significant edits ; # changed p.thres to 0.01 - Karen and Diu
diff.frequency.thres <- 0.1 # cutoff value for difference in editing frequency
selected.annotation <- c("utr3", "cds", "utr5") # genomic region

# define significance criteria for gene expression
fpkm.thres <- 5

# filter edit data based on previous significance criteria ; added DCD.fpkm and MIG.fpkm and altered diff.frequency filter - Karen and Diu
filtered.edit.lsc <- subset(raw.edit.lsc, 
                            p.adj < p.thres & 
                              annotation %in% selected.annotation &
                              diff.frequency >= diff.frequency.thres & 
                              ADA.fpkm >= fpkm.thres &
                              DCD.fpkm >= fpkm.thres &
                              MIG.fpkm >= fpkm.thres)
filtered.edit.lsk <- subset(raw.edit.lsk,
                            p.adj < p.thres & 
                              annotation %in% selected.annotation &
                              diff.frequency >= diff.frequency.thres & 
                              ADA.fpkm >= fpkm.thres &
                              DCD.fpkm >= fpkm.thres &
                              MIG.fpkm >= fpkm.thres)

filtered.gene.lsc <- subset(raw.gene.lsc,
                            ADA.fpkm > fpkm.thres)
filtered.gene.lsk <- subset(raw.gene.lsk,
                            ADA.fpkm > fpkm.thres)

# make a venn diagram for filtered data
vd <- venn.diagram(
  x = list(
    "LSK" = filtered.edit.lsk$site.id,
    "LSC" = filtered.edit.lsc$site.id
  ),
  filename = "./editSiteComparison.tiff",
  col = "transparent",
  fill = c("cornflowerblue", "green"),
  alpha = 0.50,
  label.col = c("blue", "blue", "blue"),
  scaled = TRUE,
  ext.text = TRUE,
  ext.line.lwd = 2,
  ext.dist = -0.15,
  ext.length = 0.9,
  ext.pos = -4,
  inverted = FALSE,
  cex = 1.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen"),
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  rotation.degree = 0
)

# plot distributions of edit frequency between lsc/lsk datasets
edit.site.overlap <- merge(filtered.edit.lsc[,c("site.id", "entrez.id", "diff.frequency")], 
                           filtered.edit.lsk[,c("site.id", "entrez.id", "diff.frequency")], 
                           by=c("site.id", "entrez.id"),
                           suffixes=c(".lsc", ".lsk"))
# rename columns
names(edit.site.overlap)[names(edit.site.overlap == "diff.frequency.lsc")] <- "lsc"
names(edit.site.overlap)[names(edit.site.overlap == "diff.frequency.lsk")] <- "lsk"

df <- gather(edit.site.overlap, key="cell.type", value="diff.frequency",
             diff.frequency.lsc, diff.frequency.lsk)

p1 <- ggplot(df, aes(x=cell.type, y=diff.frequency)) + 
  geom_violin(aes(fill=cell.type)) + 
  geom_boxplot(width=0.2) + 
  labs(title="HyperTRIBE Edit Frequency", x="Cell type", y="") +
  theme_light() +
  theme(legend.position="none")

# what is the average fold change at each shared editing site?
edit.site.overlap$foldchange <- log2(edit.site.overlap$diff.frequency.lsc / edit.site.overlap$diff.frequency.lsk)

p2 <- ggplot(data=edit.site.overlap, aes(edit.site.overlap$foldchange)) + 
  geom_histogram(fill="royalblue4", 
                 alpha=0.5,
                 bins=40,
                 colour="white") + 
  labs(title="Log2 Fold Change in Edit Frequency (LSC/LSK)") +
  theme_bw() +
  xlab("\nLog2 Fold Change in Edit Freq (LSC/LSK)") +
  ylab("Count\n") +
  ggtitle("Shared Sites ; padj < 0.01\n") +
  theme(plot.title = element_text(size=40)) + # make title font size 40
  theme(axis.text=element_text(size=30, color="black"), axis.title=element_text(size=30, color="black")) + # make axis labels bigger
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black", size=1)) # remove top and right plot border line

p2 # Need to enter this command to be able to plot histogram - Karen and Diu

# what is the change in edits per gene?
# get edits per gene in filtered dataset
lsc.edited.genes <- as.data.frame(base::table(filtered.edit.lsc$entrez.id))
colnames(lsc.edited.genes) <- c("entrez.id", "n.sites")
lsk.edited.genes <- as.data.frame(base::table(filtered.edit.lsk$entrez.id))
colnames(lsk.edited.genes) <- c("entrez.id", "n.sites")

# get edits in both genes
merge.edited.genes <- merge(lsc.edited.genes, lsk.edited.genes, 
                            by=c("entrez.id"),
                            suffix=c(".lsc", ".lsk"))

# add an extra column (cell.type) giving lsc/lsk, and combined n.sites into a single column
merge.edited.genes.toplot <- gather(merge.edited.genes, 
                                    key = "cell.type",
                                    value = "n.sites",
                                    n.sites.lsc,
                                    n.sites.lsk)
# rename cell types
merge.edited.genes.toplot$cell.type <- sapply(merge.edited.genes.toplot$cell.type, 
                                              function(x) stringr::str_sub(x, -3, -1))

p3 <- ggplot(merge.edited.genes.toplot, aes(x=cell.type, y=n.sites)) + 
  geom_violin(aes(fill=cell.type)) + 
  geom_boxplot(width=0.2) + 
  labs(title="HyperTRIBE Edit Frequency", x="Cell type", y="") +
  theme_light() +
  theme(legend.position="none")