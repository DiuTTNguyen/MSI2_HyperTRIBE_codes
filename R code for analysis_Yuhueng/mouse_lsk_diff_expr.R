library( DESeq2 )
library( BiocParallel )

se <- readRDS( "rnaseq_read_count_entrez_id.rds" )
colnames( se ) <- substr( colnames(se ), 1, 12 )
condition <- substr( colnames(se), 8, 10 )
colData( se )$condition <- factor( condition )
dds <- DESeqDataSet( se, design = ~condition )
fpkm.mat <- fpkm( dds, robust = T )

rowmask <- rowMeans( fpkm.mat ) > median( rowMeans( fpkm.mat ) )
colmask <- condition %in% c( 'ADA', 'DCD' )
dds1 <- dds[ rowmask, colmask ]

dds1$condition <- droplevels( dds1$condition )
dds1 <- DESeq( dds1, betaPrior = F, parallel = T, BPPARAM=MulticoreParam(4) )
res  <- results( dds1 )

colmask <- condition %in% c( 'ADA', 'MIG' )
dds2 <- dds[ rowmask, colmask ]

dds2$condition <- droplevels( dds2$condition )
dds2 <- DESeq( dds2, betaPrior = F, parallel = T, BPPARAM=MulticoreParam(4) )
res  <- results( dds2 )

pdf( 'mouse_lsk_ada_vs_dcd.pdf', 6, 5, useDingbats = F )
plotMA( dds1, main = "Gene expression: log2(DCD/ADA)" )
abline( h = c( -.5, .5 ), lty = 3 )
dev.off()
pdf( 'mouse_lsk_ada_vs_mig.pdf', 6, 5, useDingbats = F )
plotMA( dds2, main = "Gene expression: log2(MIG/ADA)" )
abline( h = c( -.5, .5 ), lty = 3 )
dev.off()
