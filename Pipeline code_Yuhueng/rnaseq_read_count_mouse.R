library(VariantAnnotation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(BiocParallel)
library(parallel)
library(GenomicAlignments)
library(Rsamtools)
library(reshape2)
library(openxlsx) 

bam.files <- Sys.glob("Sample*[0-9].dedupped.bam")

filter.vcf <- function( vcf.file, exons, genome ){
    vcf <- readVcf( vcf.file, genome )
# Filter all SNPs found in dbSNP
    vcf <- vcf[ !info(vcf)$DB ]
# Only keep A/G mutations in forward strand and T/C mutations in reverse strand
    mask1 <- rowRanges( vcf )$REF == 'A' & sapply( rowRanges( vcf )$ALT, function( seq ) return( 'G' %in% seq ) )
    mask2 <- rowRanges( vcf )$REF == 'T' & sapply( rowRanges( vcf )$ALT, function( seq ) return( 'C' %in% seq ) )
    gr1 <- rowRanges( vcf )[ mask1 ]
    gr2 <- rowRanges( vcf )[ mask2 ]
    strand(gr1) <- '+'
    strand(gr2) <- '-'
    a2g.snp <- c( gr1, gr2 )
# Only keep A/G mutations within exons
    a2g.snp <- a2g.snp[ countOverlaps( a2g.snp, exons ) > 0 ]
    return(a2g.snp)
}

vcf.files <- Sys.glob( "*FinalR.vcf" )
# Load all the exons from UCSC gene models
exbygene <- exonsBy( TxDb.Mmusculus.UCSC.mm10.knownGene, "gene" )
l <- mclapply( vcf.files, filter.vcf, exbygene, "mm10", mc.cores = length(vcf.files) )
names(l) <- basename( vcf.files )
grl <- GRangesList(l)
saveRDS( grl, file = "a2g_snp_filtered.rds" )

# grl <- readRDS("a2g_snp_filtered.rds")
a2g.snp <- sort( unique( unlist( grl, use.names = F ) ) ) 

pileup.res <- mclapply( bam.files, function( bf ){
    bf <- BamFile( bf )
    res <- pileup( bf, a2g.snp, 
        scanBamParam=ScanBamParam( flag = scanBamFlag( hasUnmappedMate=F, isProperPair=T, isDuplicate=F ), which = a2g.snp ), 
        pileupParam=PileupParam( distinguish_strands=F, min_base_quality=10, max_depth=1e4 ) )
    return( res ) }, mc.cores = length( bam.files ) )
names( pileup.res ) <- basename( bam.files )
saveRDS( pileup.res, file = "pileup_res.rds" )

# pileup.res <- lapply( pileup.res, dcast, which_label ~ nucleotide, value.var = 'count', fill = 0, drop = F )

# if( F ){
pileup.res <- lapply( pileup.res, function( df ){
    df <- dcast( df, which_label ~ nucleotide, value.var = 'count', fill = 0, drop = F )
    df$strand <- as.character( strand( a2g.snp ) )
    snp.id <- sprintf( "%s:%d-%d", seqnames( a2g.snp ), start( a2g.snp ), start( a2g.snp ) )
    stopifnot( all( snp.id == df$which_label ) )
    rownames( df ) <- df$which_label
    df <- split( df, df$strand ) 
    df$`+` <- data.frame( ref = 'A', alt = 'G', ref.count = df$`+`$A, alt.count = df$`+`$G, row.names = rownames( df$`+` ) )
    df$`-` <- data.frame( ref = 'T', alt = 'C', ref.count = df$`-`$T, alt.count = df$`-`$C, row.names = rownames( df$`-` ) )
    df <- rbind( df$`+`, df$`-` )[ snp.id, ]
    return( df )
} )

# Most simplistic allele annotation
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

# Should we go back and check the edits found in dbSNP? Maybe some of them are real after all?
# Add back dbSNP annotation
# dbsnp.gr <- endoapply( grl, function( gr ) gr[ grep( '^rs', names(gr) )] )
# dbsnp.gr <- sort( unique( unlist( dbsnp.gr, use.names = F ) ) )

# anno$dbsnp <- ""
# ol <- findOverlaps( a2g.snp.subset, dbsnp.gr )
# anno$dbsnp[ queryHits( ol ) ] <- names( dbsnp.gr )[ subjectHits( ol ) ]

mcols( a2g.snp.subset ) <- anno


pileup.res <- lapply( pileup.res, function( df, a2g.snp ){
    snp.id <- sprintf( "%s:%d-%d", seqnames( a2g.snp ), start( a2g.snp ), start( a2g.snp ) )
    df <- df[ snp.id, ]
    df <- cbind( as.data.frame( a2g.snp), df )
    df$end <- NULL
    df$width <- NULL
    return( df )
}, a2g.snp.subset )
anno.df <- pileup.res[[1]][(1:(ncol(pileup.res[[1]])-2))]
pileup.res <- lapply( pileup.res, function( df ) df[-(1:(ncol(pileup.res[[1]])-2))] )
pileup.res <- do.call( "cbind", pileup.res )
pileup.res <- cbind( anno.df, pileup.res )

write.xlsx( pileup.res, file = "snp_counts_dedupped.xlsx" )

# }
