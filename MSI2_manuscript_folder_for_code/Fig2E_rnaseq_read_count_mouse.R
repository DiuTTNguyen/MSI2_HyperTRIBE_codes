library(VariantAnnotation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(BiocParallel)
library(parallel)
library(GenomicAlignments)
library(Rsamtools)
library(reshape2)
library(openxlsx) # xlsx package completely unusable

register(MulticoreParam(workers = 6))

bam.files <- Sys.glob("*.split.bam")

ebg1 <- exonsBy( TxDb.Mmusculus.UCSC.mm10.knownGene, by="gene" )
se1 <- summarizeOverlaps(features=ebg1, reads=bam.files,
        mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE )
saveRDS( se1, file = "rnaseq_read_count_entrez_id.rds" )


filter.vcf <- function( vcf.file, exons, genome ){
    vcf <- readVcf( vcf.file, genome )
# Filter all SNPs found in dbSNP
    # vcf <- vcf[ !info(vcf)$DB ]
# Only keep A/G mutations in positive strand and T/C mutations in reverse strand
    mask1 <- rowRanges( vcf )$REF == 'A' & sapply( rowRanges( vcf )$ALT, function( seq ) return( 'G' %in% seq ) ) # obtain A to G SNPs
    mask2 <- rowRanges( vcf )$REF == 'T' & sapply( rowRanges( vcf )$ALT, function( seq ) return( 'C' %in% seq ) ) # obtain T to C SNPs
    gr1 <- rowRanges( vcf )[ mask1 ]
    gr2 <- rowRanges( vcf )[ mask2 ]
    strand(gr1) <- '+' # Assign A to G as + strand
    strand(gr2) <- '-' # Assign T to C as - strand
    a2g.snp <- c( gr1, gr2 )
# Only keep A/G mutations within exons
    a2g.snp <- a2g.snp[ countOverlaps( a2g.snp, exons ) > 0 ] # Select a2g SNPs that occur in exons.
    return(a2g.snp)
}

# if( F ){
vcf.files <- Sys.glob( "*FinalR.vcf" )
# Load all the exons from UCSC gene models
exbygene <- exonsBy( TxDb.Mmusculus.UCSC.mm10.knownGene, "gene" )
l <- mclapply( vcf.files, filter.vcf, exbygene, "mm10", mc.cores = length(vcf.files) ) # Obtain a2g, t2c snps in exons
names(l) <- basename( vcf.files )
grl <- GRangesList(l)
saveRDS( grl, file = "a2g_snp_filtered.rds" )
# }

grl <- readRDS("a2g_snp_filtered.rds")
a2g.snp <- sort( unique(  unlist( grl, use.names = F ) ) ) # Does this de-duplicate the dataframe?
pileup.res <- mclapply( bam.files, function( bf ){
    bf <- BamFile( bf )
    # pileup = Counts reads overlapping specific genomic position and assigns counts to ref(A,T) and alt(G,C) - Karen
    # pileup will also eliminate variants with zero reads. If there are 0 ref counts, the position for A is eliminated. - Karen
    res <- pileup( bf, a2g.snp, 
        scanBamParam=ScanBamParam( flag = scanBamFlag( hasUnmappedMate=F, isProperPair=T, isDuplicate=F ), which = a2g.snp ), 
        pileupParam=PileupParam( distinguish_strands=F, min_base_quality=10, max_depth=1e4 ) )
    return( res ) }, mc.cores = length( bam.files ) )
names( pileup.res ) <- basename( bam.files )
# pileup.res <- lapply( pileup.res, dcast, which_label ~ nucleotide, value.var = 'count', fill = 0, drop = F )

# if( F ){
pileup.res <- lapply( pileup.res, function( df ){
    df <- dcast( df, which_label ~ nucleotide, value.var = 'count', fill = 0, drop = F )
    df$strand <- as.character( strand( a2g.snp ) )
    # snp.id = chr and position of SNP
    snp.id <- sprintf( "%s:%d-%d", seqnames( a2g.snp ), start( a2g.snp ), start( a2g.snp ) )
    stopifnot( all( snp.id == df$which_label ) )
    rownames( df ) <- df$which_label
    df <- split( df, df$strand ) 
    df$`+` <- data.frame( ref = 'A', alt = 'G', ref.count = df$`+`$A, alt.count = df$`+`$G, row.names = rownames( df$`+` ) )
    df$`-` <- data.frame( ref = 'T', alt = 'C', ref.count = df$`-`$T, alt.count = df$`-`$C, row.names = rownames( df$`-` ) )
    df <- rbind( df$`+`, df$`-` )[ snp.id, ]
    return( df )
} )
saveRDS( pileup.res, file = "pileup_res.rds" )

# Most simplistic allele annotation
# exbygene <- exonsBy( TxDb.Hsapiens.UCSC.hg19.knownGene, "gene" )
idx <- which( countOverlaps( a2g.snp, exbygene ) == 1 )
a2g.snp.subset <- a2g.snp[idx]
ol  <- findOverlaps( a2g.snp.subset, exbygene )
entrez.id <- names( exbygene )[ subjectHits(ol) ]
gene.symbol <- mget( entrez.id, org.Mm.egSYMBOL, ifnotfound = NA )
stopifnot( all( elementNROWS(gene.symbol) == 1 ) )
gene.symbol <- unlist( gene.symbol )
anno <- DataFrame( entrez.id, gene.symbol )

anno$annotation <- "cds"
utr3bytx <- threeUTRsByTranscript(TxDb.Mmusculus.UCSC.mm10.knownGene)
utr5bytx <- fiveUTRsByTranscript(TxDb.Mmusculus.UCSC.mm10.knownGene)
anno$annotation[countOverlaps( a2g.snp.subset, utr5bytx ) > 0] <- 'utr5'
anno$annotation[countOverlaps( a2g.snp.subset, utr3bytx ) > 0] <- 'utr3'

# Add back dbSNP annotation
dbsnp.gr <- endoapply( grl, function( gr ) gr[ grep( '^rs', names(gr) )] )
dbsnp.gr <- sort( unique( unlist( dbsnp.gr, use.names = F ) ) )

anno$dbsnp <- ""
ol <- findOverlaps( a2g.snp.subset, dbsnp.gr )
anno$dbsnp[ queryHits( ol ) ] <- names( dbsnp.gr )[ subjectHits( ol ) ]

mcols( a2g.snp.subset ) <- anno


pileup.res <- lapply( pileup.res, function( df, a2g.snp ){
    snp.id <- sprintf( "%s:%d-%d", seqnames( a2g.snp ), start( a2g.snp ), start( a2g.snp ) )
    df <- df[ snp.id, ]
    df <- cbind( as.data.frame( a2g.snp), df )
    df$end <- NULL
    df$width <- NULL
    return( df )
}, a2g.snp.subset )

# saveRDS( pileup.res, file = "pileup_res.rds" )

df <- do.call( "cbind", lapply( pileup.res, function(df) df[-(1:9)] ) )
df <- cbind( pileup.res[[1]][1:9], df )
write.xlsx( df, file = "Mouse_HSPC_snps_count.xlsx" )


# }

# Karen's notes
# Seems Yuheng is taking the union of all variants called from all samples. This is why he uses pileup function to call counts from the bam files rather than subsetting the read counts from the VCF. Because if a sample doesn't have a variant that another sample has, then it's alternative read counts will be zero (and reference read count is > 0), GATK will never call it as a variant, and therefore GATK will never report the ref/alt read counts. By calling pileup function, it will produce read counts from bam file.










