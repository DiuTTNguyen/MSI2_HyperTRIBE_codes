require(GenomicFeatures)
require(ggplot2)
require(openxlsx)
require(org.Hs.eg.db)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(VennDiagram)

setwd(file.path("~", "..", "Dropbox (MIT)", "MIT", "rotation_2"))

# load transcript database
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# get entrez id -> gene name mapping
x <- org.Hs.eg.db::org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

# load genes (returns GRanges object)
hg19.genes <- GenomicFeatures::genes(txdb)

# start from GenomicRanges (code for making this on lilac cluster)
data.dir <- file.path(".", "hypertribe_data")
peak.data.fname <- file.path(data.dir, "peakdata.rds")
gr.fname <- file.path(data.dir, "CLIP_genomic_ranges.RData")

# molm13 raw dataframe (for unfiltered data)
molm13.edit.fname <- file.path(data.dir, "molm13_snp_counts_dedupped_significant.csv")

# functions
CalculateDistance <- function(query.idx, subject.idx, strand) {
  if(length(query.idx) == 0 || length(subject.idx) == 0 || length(strand) == 0) {
    return(NA)
  }
  if(as.character(strand)=="+"){
    return(query.idx - subject.idx)
  } else {
    return(subject.idx - query.idx)
  }
}

# loads the following Granges:
# k562.signif.gr
# molm13.signif.gr
# nb4.signif.gr
# peaks.signif.gr
load(gr.fname)

# peaks.data contains peaks found by CLIPanalyze
peaks.data <- readRDS(peak.data.fname)

# dataframe of molm13 HyperTRIBE hits
molm13.edit.df <- read.csv(molm13.edit.fname)

# apply new filters
ada.fpkm.thres <- 5
dcd.fpkm.thres <- 5
mig.fpkm.thres <- 5
diff.freq.thres <- 0.1

# subset molm13 dataframe based on new filters (previously only applied ada.fpkm and diff.freq filters)
molm13.filtered.df <- subset(molm13.edit.df, 
                             ADA.fpkm > ada.fpkm.thres &
                               DCD.fpkm > dcd.fpkm.thres &
                               MIG.fpkm > mig.fpkm.thres &
                               diff.frequency > diff.freq.thres)

# create Genomic Ranges object of MOLM13 data after more stringent filters
molm13.filtered.df$start <- molm13.filtered.df$pos
molm13.filtered.df$end <- molm13.filtered.df$pos
molm13.filtered.gr <- GenomicRanges::makeGRangesFromDataFrame(molm13.filtered.df)

# get subsets of K562 and NB4
nb4.raw.df <- openxlsx::read.xlsx(file.path("hypertribe_data", "nb4_clip.xlsx"), sheet = 7, startRow = 3, colNames = TRUE)
nb4.peaks.df <- rbind(nb4.raw.df[,2:8], nb4.raw.df[,9:15])
colnames(nb4.peaks.df) <- c("chr", "start", "end", "ensembl.id", "p", "strand", "gene.symbol")

# order by P value
nb4.filt.df <- nb4.peaks.df[order(nb4.peaks.df$p),]
nb4.filt.df <- nb4.filt.df[!duplicated(nb4.filt.df$gene.symbol),]
nb4.filt.df <- nb4.filt.df[order(nb4.filt.df$p),]
# nb4.small.gr <- GenomicRanges::makeGRangesFromDataFrame(nb4.filt.df[1:3000,])

# peaks small
peaks.small.gr <- peaks.data$peaks[order(peaks.data$peaks@elementMetadata$padj)]
# peaks.small.gr <- peaks.small.gr[1:3000]

#molm13 small
molm13.small.df <- molm13.filtered.df[order(molm13.filtered.df$p.adj),]
molm13.small.df <- molm13.small.df[!duplicated(molm13.small.df$entrez.id),]

# query TxDb for hits
# nb4.hits <- GenomicRanges::findOverlaps(nb4.signif.gr, hg19.genes)
# k562.hits <- GenomicRanges::findOverlaps(peaks.data$peaks, hg19.genes)
# as per Diu's suggestion, only take top 3000 genes
nb4.hits <- GenomicRanges::findOverlaps(nb4.small.gr, hg19.genes)
k562.hits <- GenomicRanges::findOverlaps(peaks.small.gr, hg19.genes)
molm13.hits <- GenomicRanges::findOverlaps(molm13.filtered.gr, hg19.genes)

# generate dataframe of Entrez ID's and gene symbols
nb4.hits.gr <- hg19.genes[subjectHits(nb4.hits)]
k562.hits.gr <- hg19.genes[subjectHits(k562.hits)]
molm13.hits.gr <- hg19.genes[subjectHits(molm13.hits)]

# add new column with gene symbol as element metadata
nb4.hits.gr@elementMetadata$gene_symbol <- xx[nb4.hits.gr@elementMetadata$gene_id]
k562.hits.gr@elementMetadata$gene_symbol <- xx[k562.hits.gr@elementMetadata$gene_id]
molm13.hits.gr@elementMetadata$gene_symbol <- xx[molm13.hits.gr@elementMetadata$gene_id]

# deduplicate
nb4.hits.gr <- nb4.hits.gr[!duplicated(nb4.hits.gr)]
k562.hits.gr <- k562.hits.gr[!duplicated(k562.hits.gr)]
molm13.hits.gr <- molm13.hits.gr[!duplicated(molm13.hits.gr)]

# make dataframes
nb4.hits.df <- GenomicRanges::as.data.frame(nb4.hits.gr)
k562.hits.df <- GenomicRanges::as.data.frame(k562.hits.gr)
molm13.hits.df <- GenomicRanges::as.data.frame(molm13.hits.gr)

# write some CSV files
write.csv(nb4.hits.df, file = file.path("clip", "nb4_genes.csv"))
write.csv(k562.hits.df, file = file.path("clip", "k562_genes.csv"))
write.csv(molm13.hits.df, file = file.path("clip", "molm13_genes.csv"))

# new filtered df
write.csv(nb4.hits.df, file = file.path("clip", "nb4_short_genes.csv"))
write.csv(k562.hits.df, file = file.path("clip", "k562_short_genes.csv"))
write.csv(molm13.small.df, file = file.path("clip", "molm13_short_genes.csv"))

# NOTE CHANGED MOLM13 GENE IDS
nb4.molm13.vd <- VennDiagram::venn.diagram(
  x = list(
    "NB4" = nb4.hits.gr$gene_id,
    "MOLM13" = molm13.small.df$entrez.id
  ),
  filename = file.path("plots", "nb4_molm13_venn.tiff"),
  col = "transparent",
  fill = c("cornflowerblue", "green"),
  alpha = 0.50,
  label.col = rep("blue", 3),
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

k562.molm13.vd <- VennDiagram::venn.diagram(
  x = list(
    "K562" = k562.hits.gr$gene_id,
    "MOLM13" = molm13.small.df$entrez.id
  ),
  filename = file.path("plots", "k562_molm13_venn.tiff"),
  col = "transparent",
  fill = c("cornflowerblue", "green"),
  alpha = 0.50,
  label.col = rep("blue", 3),
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

# get indices of nearest neighbors
molm13.nb4.nearest <- GenomicRanges::nearest(molm13.signif.gr, nb4.signif.gr)
molm13.k562.nearest <- GenomicRanges::nearest(molm13.signif.gr, k562.signif.gr)
molm13.peaks.nearest <- GenomicRanges::nearest(molm13.signif.gr, peaks.signif.gr)
molm13.peaks.unfiltered.nearest <- GenomicRanges::nearest(molm13.signif.gr, peaks.data$peaks)

# get indices of nearest neighbors (the same shit for more stringently filtered data)
molm13.filtered.nb4.nearest <- GenomicRanges::nearest(molm13.filtered.gr, nb4.signif.gr)
molm13.filtered.k562.nearest <- GenomicRanges::nearest(molm13.filtered.gr, k562.signif.gr)
molm13.filtered.peaks.nearest <- GenomicRanges::nearest(molm13.filtered.gr, peaks.data$peaks)

# replace all NA values with 0
molm13.nb4.nearest[is.na(molm13.nb4.nearest)] <- 0
molm13.k562.nearest[is.na(molm13.k562.nearest)] <- 0
molm13.peaks.nearest[is.na(molm13.peaks.nearest)] <- 0

# reaplce all NA values with 0 (same shit for more stringently filtered data)
molm13.filtered.nb4.nearest[is.na(molm13.filtered.nb4.nearest)] <- 0
molm13.filtered.k562.nearest[is.na(molm13.filtered.k562.nearest)] <- 0
molm13.filtered.peaks.nearest[is.na(molm13.filtered.peaks.nearest)] <- 0

# get NB4/K562 distances (publicly available datasets)
molm13.nb4.distance <- mapply(function(q.idx, s.idx, strand) CalculateDistance(q.idx, s.idx, strand),
                              as.list(molm13.signif.gr@ranges@start),
                              as.list(nb4.signif.gr[c(molm13.nb4.nearest)]@ranges@start),
                              as.list(molm13.signif.gr@strand))
# Commented out because using new CLIPanalyze pipeline
# molm13.k562.distance <- mapply(function(q.idx, s.idx, strand) CalculateDistance(q.idx, s.idx, strand),
#                               as.list(molm13.signif.gr@ranges@start),
#                               as.list(k562.signif.gr[c(molm13.k562.nearest)]@ranges@start),
#                               as.list(molm13.signif.gr@strand))

# K562 peaks distance (from CLIPanalyze)
molm13.peaks.distance <- mapply(function(q.idx, s.idx, strand) CalculateDistance(q.idx, s.idx, strand),
                                as.list(molm13.signif.gr@ranges@start),
                                as.list(peaks.data$peaks[c(molm13.peaks.unfiltered.nearest)]@ranges@start),
                                as.list(molm13.signif.gr@strand))

# again for more stringent filters (LOL this code is horrible.)
# get NB4/K562 distances (publicly available datasets)
molm13.filtered.nb4.distance <- mapply(function(q.idx, s.idx, strand) CalculateDistance(q.idx, s.idx, strand),
                              as.list(molm13.filtered.gr@ranges@start),
                              as.list(nb4.signif.gr[c(molm13.filtered.nb4.nearest)]@ranges@start),
                              as.list(molm13.filtered.gr@strand))
# Commented out because using new CLIPanalyze pipeline
# molm13.filtered.k562.distance <- mapply(function(q.idx, s.idx, strand) CalculateDistance(q.idx, s.idx, strand),
#                                as.list(molm13.filtered.gr@ranges@start),
#                                as.list(k562.signif.gr[c(molm13.filtered.k562.nearest)]@ranges@start),
#                                as.list(molm13.filtered.gr@strand))

# K562 peaks distance (from CLIPanalyze)
molm13.filtered.peaks.distance <- mapply(function(q.idx, s.idx, strand) CalculateDistance(q.idx, s.idx, strand),
                                as.list(molm13.filtered.gr@ranges@start),
                                as.list(peaks.data$peaks[c(molm13.filtered.peaks.nearest)]@ranges@start),
                                as.list(molm13.filtered.gr@strand))

# analysis of K562 peaks vs NB4 peaks (K562 peaks called with CLIPanalyze)
k562.nb4.dist <- GenomicRanges::distanceToNearest(peaks.data$peaks, nb4.signif.gr)
k562.nb4.df <- data.frame(chr = peaks.data$peaks@seqnames,
                          pos = peaks.data$peaks@ranges@start,
                          dist = k562.nb4.dist@elementMetadata$distance)
k562.nb4.hist <- ggplot(data=k562.nb4.df, aes(k562.nb4.df$dist)) + 
  geom_histogram(breaks=seq(0,4000, by=20),
                 col="white", 
                 fill="blue", 
                 alpha=0.2) + 
  labs(title="Distance from K562 peak to nearest NB4 CLIP site") +
  labs(x="Distance (bp)") +
  labs(y="Count")


# create df with all of the distances
distance.df <- data.frame(chr <- molm13.signif.gr@seqnames,
                          strand <- molm13.signif.gr@strand,
                          pos <- molm13.signif.gr@ranges@start,
                          nb4.dist <- molm13.nb4.distance,
                          # k562.dist <- molm13.k562.distance,
                          peaks.dist <- molm13.peaks.distance)

# create df with all of the distances (filtered MOLM13)
filtered.distance.df <- data.frame(chr <- molm13.filtered.gr@seqnames,
                          strand <- molm13.filtered.gr@strand,
                          pos <- molm13.filtered.gr@ranges@start,
                          nb4.dist <- molm13.filtered.nb4.distance,
                          # k562.dist <- molm13.filtered.k562.distance,
                          peaks.dist <- molm13.filtered.peaks.distance)

# NB4histogram
nb4.hist <- ggplot(data=distance.df, aes(distance.df$nb4.dist)) + 
  geom_histogram(breaks=seq(-1000,1000, by=50), 
                 col="white", 
                 fill="blue", 
                 alpha=0.2) + 
  labs(title="Distance to nearest NB4 CLIP site") +
  labs(x="Distance (bp)") +
  labs(y="Count")

# K562 histogram (using old CLIP hits) (commented out because using new CLIPanalyze pipeline)
# k562.hist <- ggplot(data=distance.df, aes(distance.df$k562.dist)) + 
#   geom_histogram(breaks=seq(-100000, 100000, by=5000),
#                  col="white", 
#                  fill="blue", 
#                  alpha=0.2) + 
#   labs(title="Distance to nearest NB4 CLIP site") +
#   labs(x="Distance (bp)") +
#   labs(y="Count")

# K562 histogram (using new CLIPanalyze algo)
peaks.hist <- ggplot(data=distance.df, aes(distance.df$peaks.dist)) + 
  geom_histogram(breaks=seq(-1000,1000, by=50),
                 col="white", 
                 fill="blue", 
                 alpha=0.2) + 
  labs(title="Distance to nearest K562 CLIP site") +
  labs(x="Distance (bp)") +
  labs(y="Count")

# same histograms, filtered data
# NB4histogram
nb4.filtered.hist <- ggplot(data=filtered.distance.df, aes(filtered.distance.df$nb4.dist)) + 
  geom_histogram(breaks=seq(-1000,1000, by=50), 
                 col="white", 
                 fill="blue", 
                 alpha=0.2) + 
  labs(title="Distance to nearest NB4 CLIP site") +
  labs(x="Distance (bp)") +
  labs(y="Count")

# K562 histogram (using old CLIP hits) (commented out, using new CLIPanalyze pipeline)
# k562.filtered.hist <- ggplot(data=filtered.distance.df, aes(filtered.distance.df$k562.dist)) + 
#   geom_histogram(breaks=seq(-100000, 100000, by=5000),
#                  col="white", 
#                  fill="blue", 
#                  alpha=0.2) + 
#   labs(title="Distance to nearest NB4 CLIP site") +
#   labs(x="Distance (bp)") +
#   labs(y="Count")

# K562 histogram (using new CLIPanalyze algo)
peaks.filtered.hist <- ggplot(data=filtered.distance.df, aes(filtered.distance.df$peaks.dist)) + 
  geom_histogram(breaks=seq(-1000,1000, by=50),
                 col="white", 
                 fill="blue", 
                 alpha=0.2) + 
  labs(title="Distance to nearest K562 CLIP site") +
  labs(x="Distance (bp)") +
  labs(y="Count")
