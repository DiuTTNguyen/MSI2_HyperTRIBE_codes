# Process HSC data separately 
library(openxlsx)
library(parallel)
library(stringr)

library(bbmle)
library(VGAM)

library(DESeq2)

compartments <- c( "MPP2", "MPP4", "LT", "ST" )
df.results <- list()

for( i in compartments ){
df <- read.xlsx( "mouse_hsc_mettl3_snp_counts_dedupped_full.xlsx", sheet = i )

ref.counts <- df[,grep('ref.count', colnames(df))]
alt.counts <- df[,grep('alt.count', colnames(df))]
alt.freq <- alt.counts / ( ref.counts + alt.counts )
df.stats <- data.frame( 
                       diff.frequency = rowMeans(alt.freq[,1:3],na.rm=T) - rowMeans(alt.freq[,-(1:3)],na.rm=T),
                       ADA.frequency = rowMeans( alt.freq[,1:3],na.rm=T ),
                       DCD.frequency = rowMeans( alt.freq[,4:6],na.rm=T ),
                       MIG.frequency = rowMeans( data.frame( alt.freq[,-(1:6)]),na.rm=T ) )

rowmask <- rowSums( alt.counts != 0 ) > 1 &
( rowMeans( ref.counts[1:3] ) + rowMeans( alt.counts[1:3] ) >= 5 ) &
( rowMeans( ref.counts ) + rowMeans( alt.counts ) >= 5 ) &
!is.na(df$gene.symbol) 
df <- df[ rowmask, ]

lrt.res <- mclapply( 1:nrow(df), function( i ){ 
       counts <- as.integer( df[i,-(1:8)] )
       ref <- counts[ seq( 1, length( counts ), 2 ) ]
       alt <- counts[ seq( 2, length( counts ), 2 ) ]
       fit0 <- mle.custom.h0( ref, alt )
       fit1 <- mle.custom.h1( ref[1:3], alt[1:3], ref[-(1:3)], alt[-(1:3)] )
       p.value <- pchisq( 2 *( fit0@min - fit1@min ), 1, lower.tail=F )
       value <- fit1@min
       code <- fit1@details$convergence
       return( c( p.value = p.value, value = value, code = code ) )
}, mc.cores = 20 )
lrt.res <- data.frame( do.call( 'rbind', lrt.res ) )

df.stats <- cbind( df.stats[rowmask,], lrt.res )
rowmask <- df.stats$value > 0
df.stats <- df.stats[rowmask,]
df <- df[ rowmask, ]
df.stats$p.adj <- p.adjust( df.stats$p.value, 'BH' )
df.stats$value <- NULL
df.stats$code <- NULL
df1 <- cbind( df[ , 1:8 ], df.stats, df[ , -(1:8) ] )
df.results[[i]] <- df1
}
write.xlsx( df.results, file = "mouse_hsc_snp_counts_dedupped_new.xlsx" )
# Use FDR < 10% cutoff 
df.results <- list()
for( i in compartments ){
df <- read.xlsx( "mouse_hsc_snp_counts_dedupped_new.xlsx", sheet = i )
df1 <- df[ which( df$p.adj < .1 & df$diff.frequency > 0 ), ]
gene.num.edits <- table( df1$entrez.id )
se <- readRDS( "rnaseq_read_count_entrez_id.rds" )
se <- se[ ,colData(se)$cell == i ]
fpkm.mat <- fpkm( DESeqDataSet( se, ~condition ) )
fpkm.mat <- data.frame( ADA.fpkm = rowMeans( fpkm.mat[,1:2]),
                       DCD.fpkm = rowMeans( fpkm.mat[,3:4]),
                       MIG.fpkm = rowMeans( data.frame(fpkm.mat[,-(1:4)]) ) )
fpkm.mat <- data.frame( entrez.id = names(gene.num.edits), gene.num.edits = as.integer(gene.num.edits), fpkm.mat[ names( gene.num.edits ), ] )
df1 <- merge( df1, fpkm.mat, by = 'entrez.id' )
df.results[[i]] <- df1
}
write.xlsx( df.results, file = "mouse_hsc_snp_counts_dedupped_new.xlsx" )
