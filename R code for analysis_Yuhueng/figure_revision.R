library(openxlsx)
library(DESeq2)
library(ggplot2)

# setwd( "~/projects/Msi2_ADAR/" )
# MA plots for four HSC compartments
# CSV files are under google drive folder DESeq2 results => HSC subpopulations
compartments <- c( 'lt', 'st', 'mpp2', 'mpp4' )
for( x in compartments ){
    df <- read.csv( sprintf( 'hsc_%sdcd_vs_%smig.csv', x, x ), row.names = 1  )
    df <- df[ order( df$padj, decreasing = T ), ]
    rowmask1 <- df$log2FoldChange > 6
    rowmask2 <- df$log2FoldChange < -6
    rowmask3 <- !( rowmask1 | rowmask2 )
    df$shape <- 'mid'
    df$shape[ rowmask1 ] <- 'hi'
    df$shape[ rowmask2 ] <- 'lo'
    df$log2FoldChange[ rowmask1 ] <- 6
    df$log2FoldChange[ rowmask2 ] <- -6

    pdf( sprintf( 'hsc_%sdcd_vs_%smig.pdf', x, x ), 5, 4, useDingbats = F )
    qplot( baseMean, log2FoldChange, data = df, geom = 'point', color = padj > .05, shape = shape ) +
        geom_hline( yintercept = 0, linetype = 3 ) + 
        scale_x_log10() + 
        scale_color_manual( values = c( 'red3', 'grey32' ) ) + 
        scale_shape_manual( values = c( 2, 20, 6 ), limits = c( 'hi', 'mid', 'lo' ) ) + 
        ylim( -6, 6 ) + 
        xlab( "mean of normalized counts" ) + ylab( "log2 fold change" ) + 
        theme_classic() + theme( legend.position = 'none' )
    dev.off()
}  


# ADA.fpkm vs ADA.frequency
df <- read.xlsx( 'molm13_snp_counts_dedupped_significant.xlsx', sheet = 1 )
pdf( "supp_fig_1_alternative.pdf", 5.5, 4, useDingbats = F )
qplot( x = ADA.frequency, y = ADA.fpkm, data = df, color = '' ) +
    geom_vline( xintercept = .1, linetype = 3 ) +
    geom_hline( yintercept = 5, linetype = 3 ) +
    scale_color_manual( values = "#336A98" ) + scale_y_log10() +
    theme_classic() + theme( legend.position = 'none' )
dev.off()

# MA plot in Supp. Fig. 1
se <- SummarizedExperiment( assays = as.matrix( df[ c( grep( 'MIG_', names(df) ), grep( 'DCD_', names(df) ) ) ] ), 
                           colData = DataFrame( condition = c( rep( 'MIG', 3 ), rep( 'DCD', 3 ) ) ) )

rowmask <- df$MIG.fpkm >= 5 & df$DCD.fpkm >= 5
dds <- DESeqDataSet( se[ rowmask, ], ~condition )
dds <- DESeq( dds, betaPrior = F )
res <- results( dds )

df <- as.data.frame( res )
df <- df[ order( df$padj, decreasing = T ), ]
rowmask1 <- df$log2FoldChange > 1
rowmask2 <- df$log2FoldChange < -1
rowmask3 <- !( rowmask1 | rowmask2 )
df$shape <- 'mid'
df$shape[ rowmask1 ] <- 'hi'
df$shape[ rowmask2 ] <- 'lo'
df$log2FoldChange[ rowmask1 ] <- 1
df$log2FoldChange[ rowmask2 ] <- -1

pdf( "supp_fig_1_ma_plot.pdf", 5, 4, useDingbats = F )
qplot( baseMean, log2FoldChange, data = df, geom = 'point', color = padj > .05, shape = shape ) +
    geom_hline( yintercept = 0, linetype = 3 ) + 
    scale_x_log10() + 
    scale_color_manual( values = c( 'red3', 'grey32' ) ) + 
    scale_shape_manual( values = c( 2, 20, 6 ), limits = c( 'hi', 'mid', 'lo' ) ) + 
    ylim( -.5, .5 ) + 
    xlab( "mean of normalized counts" ) + ylab( "log2 fold change" ) + 
    theme_classic() + theme( legend.position = 'none' )
dev.off()
