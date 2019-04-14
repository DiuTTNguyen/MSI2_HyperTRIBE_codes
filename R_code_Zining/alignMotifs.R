# check:
# motif file
# csv file
# genome
# txdb

# libraries
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Set directory
setwd("C:/Users/Zi-Ning/Dropbox (MIT)/MIT/rotation_2")

# Get transcript db
# txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# get GRanges object containing 3'UTRs (or exons)
utr.gr <- GenomicFeatures::exonsBy(txdb, by=("tx"), use.names = TRUE)
# utr.gr <- threeUTRsByTranscript(txdb, use.names = TRUE)
# utr.gr <- threeUTRsByTranscript(txdb, columns=c("exon_id"))

# CSV files containing HyperTRIBE edit data
data.path <- "hypertribe_data"
csv.fname <- "mouse_lsc_snp_counts_dedupped_significant.csv"
# csv.fname <- "mouse_lsk_snp_counts_dedupped_significant.csv"
# csv.fname <- "mouse_hsc_snp_counts_dedupped_significant_mpp2.csv"
# csv.fname <- "mouse_hsc_snp_counts_dedupped_significant_mpp4.csv"
# csv.fname <- "mouse_hsc_snp_counts_dedupped_significant_st.csv"
# csv.fname <- "mouse_hsc_snp_counts_dedupped_significant_lt.csv"
# csv.fname <- "molm13_snp_counts_dedupped_significant.csv"

edit.fname <- file.path(getwd(), data.path, csv.fname)


# motif filename
motif.dir <- file.path("motif_enrichment", "lsc")
motif.fname <- file.path(getwd(), motif.dir, "denovo", "homerResults", "motif1.motif")
# motif.dir <- file.path("motif_enrichment")
# motif.fname <- file.path(motif.dir, "msi1.motif")

# output filename
output.fname <- "lsk_with_distance.csv"

# Read csv
raw <- read.csv(edit.fname)

# settings
p.thres <- 1 # cutoff p-value for significant edits
selected.annotation <- c("cds", "utr5", "utr3") # genomic region
# genome <- BSgenome.Hsapiens.UCSC.hg19
genome <- BSgenome.Mmusculus.UCSC.mm10

# Filter relevant lines
filtered <- subset(raw, p.adj <= p.thres & annotation %in% selected.annotation)

# create gr object of filtered edit sites
gr <- GRanges(seqnames = filtered$chr, 
              ranges = IRanges(start=filtered$pos, end=filtered$pos, names=filtered$pos), 
              strand = filtered$strand)

# find UTRs that contain edit sites (hits is a hits object)
utr.gr <- unlist(utr.gr)
hits <- GenomicRanges::findOverlaps(gr, unlist(utr.gr))

# granges object of edit sites that could be mapped to a UTR
query.hits <- gr[queryHits(hits),]

# granges object of UTRs that were edited
subject.hits <- utr.gr[subjectHits(hits),]

query.hits@elementMetadata$exon_id <- subject.hits@elementMetadata$exon_id

# deduplicate query and subject hits
query.hits.deduped <- query.hits[!duplicated(query.hits)]
subject.hits.deduped <- subject.hits[!duplicated(subject.hits)]

# get actual sequences of UTRs containing edits
subject.hits.seqs <- getSeq(genome, subject.hits.deduped)

# read in motif file as MATRIX (.motif files are tab-delimited)
motif.matrix <- t(as.matrix(read.table(motif.fname, skip=1, header=FALSE, sep="\t")))

# sketchy convert frequency to count (make input PFM)
tmp <- lapply(colSums(motif.matrix),function(x) rep(x, nrow(motif.matrix)))
scale.matrix <- matrix(unlist(tmp), nrow = nrow(motif.matrix), ncol = ncol(motif.matrix))
motif.matrix <- 1000000000 * motif.matrix / scale.matrix
# tmp <- rep(1, ncol(motif.matrix))
# motif.matrix[1,] <- motif.matrix[1,] + (tmp-colSums(motif.matrix))
# motif.matrix <- motif.matrix * 1000
rownames(motif.matrix) <- c("A", "C", "G", "T")
motif.matrix <- round(motif.matrix)

# PWM constructor takes INTEGER matrix as input
storage.mode(motif.matrix) <- "integer"

# get nucleotide frequencies (need this for priors)
oligo.freq.df <- oligonucleotideFrequency(subject.hits.seqs, 1)
priors <- colSums(oligo.freq.df)/sum(sum(oligo.freq.df))

# create PWM
pwm <- PWM(motif.matrix, type="log2probratio", prior.params=priors)
matched <- lapply(subject.hits.seqs, function(x) matchPWM(pwm, x, min.score="90%", with.score=TRUE))

# get matched ranges as IRangesList object
matched.irangeslist <- IRangesList(lapply(matched, function(x) as(x, "IRanges")))

# make IRangesList
GetGenomicRanges <- function(ir, gr) {
  shifted <- shift(ir, start(ranges(gr)))
  l<-length(ir)
  gr.seqname <- rep(as.character(seqnames(gr)), l)
  gr.strand <- rep(as.character(strand(gr)), l)
  gr.exon_id <- rep(gr@elementMetadata$exon_id, l)
  if( as.character(strand(gr))=="+" ) {
    return(GRanges(seqname=gr.seqname,
                   ranges=shifted,
                   strand=gr.strand,
                   exon_id=gr.exon_id))
  } else if( as.character(strand(gr))=="-" ) {
    return(GRanges(seqname=gr.seqname,
                   ranges=reflect(shifted, ranges(gr)),
                   strand=gr.strand,
                   exon_id=gr.exon_id))
  } else {
    return(NA)
  }
}

# get subject.hits.deduped as list (I AM A N00B at r)
subject.hits.gr <- as(subject.hits.deduped, "GRangesList")

# convert matched IRanges to GRanges (due to different indexing for IRanges returned by PWM matching)
test <- mapply(GetGenomicRanges, ir=matched.irangeslist, gr=subject.hits.gr)
test.grange <- unlist(GRangesList(test))

# convert genomic ranges into data frames
test.df <- as.data.frame(test.grange, row.names = NULL)
query.df <- as.data.frame(query.hits.deduped, row.names = NULL)

# get list of all start points
# l <- sapply(levels(factor(test.df$exon_id)), function(x) test.df[test.df$exon_id==x, "start"])

GetMinimumDistance <- function(exon.id, start, end, strand) {
  q <- test.df[test.df$exon_id %in% list(exon.id), c("start", "end")]
  start.vector <- q$start
  end.vector <- q$end
  if(length(start.vector) == 0 || length(end.vector) == 0) {
    return(NA)
  }
  vector.idx <- which.min(abs(start.vector - start))
  if(strand=="+") {
    return(start - start.vector[vector.idx])
  } else if(strand=="-") {
    return(end.vector[vector.idx] - start)
  } else {
    return(NA)
  }
}

subject.hit.df <- as.data.frame(subject.hits.deduped, row.names = NULL)

# get closest start point on chromosome
shortest.dist <- mapply(function(exon.id, start, end, strand) GetMinimumDistance(exon.id, start, end, strand),
                        query.df$exon_id,
                        query.df$start,
                        query.df$end,
                        query.df$strand)

# add shortest distance to query df
query.df$dist.to.motif <- unlist(shortest.dist)

# merge with raw dataframe
dist.df <- data.frame(chr=query.df$seqnames, 
                       pos=query.df$start,
                       dist.to.motif=query.df$dist.to.motif)
raw.merge <- merge(raw, dist.df, 
                   by=c("chr", "pos"), 
                   all.x=TRUE, all.y=FALSE)

# histogram
df <- data.frame(dist=query.df$dist.to.motif)
p2 <- ggplot(data=df, aes(df$dist)) + 
  geom_histogram(breaks=seq(-500,500, by=20), 
                 col="white", 
                 fill="blue", 
                 alpha=0.2) + 
  labs(title="Distance to nearest motif") +
  labs(x="Distance (bp)") +
  labs(y="Count")

# write.csv(raw.merge, output.fname)