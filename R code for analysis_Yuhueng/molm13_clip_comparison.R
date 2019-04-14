library(ggplot2)
library(data.table)
library(Gviz)
library(openxlsx)

# TRIBE edits 
df <- read.xlsx( 'molm13_snp_counts_dedupped_new.xlsx', sheet = 1 )
tribe.edits <- GRanges( df$chr, IRanges( df$pos, df$pos ), df$strand )
mcols(tribe.edits) <- DataFrame( df[,c(1,5,6)] )

# CLIP peaks
# rep1.peaks <- import.bed("GSM1704212_MSI2_ACAGTG_ACAGTG_L008_R1.polyATrim.adapterTrim.rmDup.sorted.peaks.bed.gz")
# rep2.peaks <- import.bed("GSM1704213_MSI2_TTAGGC_TTAGGC_L008_R1.polyATrim.adapterTrim.rmDup.sorted.peaks.bed.gz")
df <- read.xlsx( 'nature17665-s2.xlsx', sheet = 7, startRow = 3 )
rep1.peaks <- GRanges( df[,2], IRanges(df[,3], df[,4]), df[,7], gene.symbol = df[,8])
idx <- which( !is.na(df[,9]) )
rep2.peaks <- GRanges( df[idx,9], IRanges(df[idx,10], df[idx,11]), df[idx,14], gene.symbol = df[idx,15] )

# File names
rep1.pos.bw <- 'GSM1704212_MSI2_ACAGTG_ACAGTG_L008_R1.polyATrim.adapterTrim.rmDup.sorted.norm.pos.bw'
rep1.neg.bw <- 'GSM1704212_MSI2_ACAGTG_ACAGTG_L008_R1.polyATrim.adapterTrim.rmDup.sorted.norm.neg.bw'
rep2.pos.bw <- 'GSM1704213_MSI2_TTAGGC_TTAGGC_L008_R1.polyATrim.adapterTrim.rmDup.sorted.norm.pos.bw'
rep2.neg.bw <- 'GSM1704213_MSI2_TTAGGC_TTAGGC_L008_R1.polyATrim.adapterTrim.rmDup.sorted.norm.neg.bw'

# Overlap between CLIP peaks & TRIBE edits
ol1 <- findOverlaps( rep1.peaks, tribe.edits )
ol2 <- findOverlaps( rep2.peaks, tribe.edits )



# Visualization
dTrack1 <- DataTrack( range = rep1.neg.bw, genome = 'hg19', type = 'l', name = 'Replicate 1' )
dTrack2 <- DataTrack( range = rep2.neg.bw, genome = 'hg19', type = 'l', name = 'Replicate 2' )
aTrack <- AnnotationTrack( range = tribe.edits, genome = 'hg19' )
plotTracks( list( aTrack, dTrack1, dTrack2 ), from =83844988, to =83846747, chromosome = 'chr16' )
chr1:26210556-26212458

extend.size <- 1000
selected.peaks <- rep1.peaks[ which( rep1.peaks$gene.symbol %in% tribe.edits$gene.symbol & width( rep1.peaks ) >= 10 ) ]
selected.peaks <- resize( selected.peaks, 2*extend.size+1, fix = 'center' ) 
selected.peaks <- selected.peaks[ countOverlaps( selected.peaks, tribe.edits ) > 0 ]
# Process +/- strands separately
selected.peaks.pos <- selected.peaks[ strand( selected.peaks ) == '+' ]
selected.peaks.neg <- selected.peaks[ strand( selected.peaks ) == '-' ]
selected.peaks <- c( selected.peaks.pos, selected.peaks.neg )

cvg.pos <- sapply( selected.peaks.pos, function( gr ){
return( unlist( import( rep1.pos.bw, format = 'bw', selection = BigWigSelection( gr ), as = "NumericList" ) ) )
} )
cvg.neg <- sapply( selected.peaks.neg, function( gr ){
return( rev( unlist( import( rep1.neg.bw, format = 'bw', selection = BigWigSelection( gr ), as = "NumericList" ) ) ) )
} )
cvg <- cbind( cvg.pos, cvg.neg )

selected.peaks.pos.ol <- findOverlaps( selected.peaks.pos, tribe.edits )
selected.peaks.neg.ol <- findOverlaps( selected.peaks.neg, tribe.edits )
edit.peak.pos.1 <- start(tribe.edits[ subjectHits(selected.peaks.pos.ol) ]) - start(selected.peaks.pos[ queryHits(selected.peaks.pos.ol) ] ) 
edit.peak.pos.2 <- end(selected.peaks.neg[ queryHits(selected.peaks.neg.ol) ] ) - end(tribe.edits[ subjectHits(selected.peaks.neg.ol) ]) 
 
edit.peak.pos <- c(edit.peak.pos.1, edit.peak.pos.2) - extend.size
dt <- data.table( peak.id = c(queryHits(selected.peaks.pos.ol), queryHits(selected.peaks.neg.ol) + max(queryHits(selected.peaks.pos.ol)) ), 
                 edit.id = c( subjectHits(selected.peaks.pos.ol), subjectHits(selected.peaks.neg.ol) ),
                 edit.peak.pos )
dt <- dt[,.(edit.peak.pos = edit.peak.pos[which(abs(edit.peak.pos) == min(abs(edit.peak.pos)) )[1]]),by=edit.id]
df <- data.frame( pos = -extend.size:extend.size, avgCov = rowMeans(cvg) )

# Try having two y axis in ggplot
ggplot() + 
geom_histogram( data = dt, aes( x = edit.peak.pos ), binwidth = 20, fill = NA, color = 'black' )  +
geom_line( data = df, aes( x = pos, y = avgCov*40 ), size = 1, color = '#f8766d' ) + 
scale_y_continuous(sec.axis = sec_axis(~./40, name = "Average CLIP coverage"))+ 
theme_classic() + xlab("Position") + ylab("TRIBE edits")
