#libraries
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)
library(IRanges)

# Set directory
setwd("C:/Users/Zi-Ning/Dropbox (MIT)/MIT/rotation_2")

# input filename
# csv.fname <- "mouse_lsc_snp_counts_dedupped_significant.csv"
# csv.fname <- "mouse_hsc_snp_counts_dedupped_significant_mpp2.csv"
# csv.fname <- "mouse_hsc_snp_counts_dedupped_significant_mpp4.csv"
# csv.fname <- "mouse_hsc_snp_counts_dedupped_significant_st.csv"
# csv.fname <- "mouse_hsc_snp_counts_dedupped_significant_lt.csv"
# csv.fname <- "mouse_lsk_snp_counts_dedupped_significant.csv"
csv.fname <- "molm13_snp_counts_dedupped_significant.csv"


# output filenames
# FORMAT: animal (mm/hs) cell (lsc/lsk/hsc) region (utr) dataset (pos/neg) 
fname.pos <- "hs_mol_utr_pos_100.fa"
fname.neg <- "hs_mol_utr_neg_100.fa"

# Read csv
raw <- read.csv(csv.fname)

# settings
p.thres <- 0.001 # cutoff p-value for significant edits
selelected.annotation <- c("utr3") # genomic region
resized.length <- 201 # size of region surrounding edited point
genome <- BSgenome.Hsapiens.UCSC.hg19

# Filter relevant lines
filtered <- subset(raw, p.adj <= p.thres & annotation == selelected.annotation)

# create gr object
gr <- GRanges(seqnames = filtered$chr, 
              ranges = IRanges(start=filtered$pos, end=filtered$pos, names=filtered$pos), 
              strand = filtered$strand)

# resize genomic ranges
resized <- resize(gr, resized.length, "center")

# reduce overlapping ranges
reduced <- reduce(resized)

# get sequences
seqs <- getSeq(genome, reduced)

# name sequences (important for making FASTA?)
names(seqs) <- paste0("pos", 1:length(seqs))

# find width distribution
seq.attr <- attributes(seqs)
range.attr <- attributes(seq.attr$ranges)
# hist(range.attr$width)

# flanking sequences
flank.width <- median(range.attr$width)
flank.beg <- flank(reduced, width = flank.width, start=TRUE, both=FALSE)
flank.end <- flank(reduced, width = flank.width, start=FALSE, both=FALSE)
flank.full <- c(flank.beg, flank.end)

# remove overlapping sequences
overlap.boolean <- countOverlaps(flank.full, reduced)
flank.nooverlaps <- flank.full[!overlap.boolean]

# get flanking sequences
flank.seqs <- getSeq(genome, flank.nooverlaps)

# name sequences (important for making FASTA?)
names(flank.seqs) <- paste0("neg", 1:length(flank.seqs))

# write sequences
writeXStringSet(seqs, fname.pos, format="fasta", append=FALSE)
writeXStringSet(flank.seqs, fname.neg, format="fasta", append=FALSE)