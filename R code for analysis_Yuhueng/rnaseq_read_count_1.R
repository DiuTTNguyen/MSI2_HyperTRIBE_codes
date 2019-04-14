library(GenomicAlignments)
library(Rsamtools)
library(BiocParallel)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

bam.files <- Sys.glob("Sample*[0-9].dedupped.bam")
ebg1 <- exonsBy( TxDb.Hsapiens.UCSC.hg19.knownGene, by="gene" )
se1 <- summarizeOverlaps(features=ebg1, reads=bam.files,
        mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE )
saveRDS( se1, file = "rnaseq_read_count_entrez_id.rds" )
