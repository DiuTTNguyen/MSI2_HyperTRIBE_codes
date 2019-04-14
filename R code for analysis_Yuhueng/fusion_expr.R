library(GenomicAlignments)
library(stringr)
library(ggplot2)

setwd('~/projects/Mouse_LSC')
l <- readRDS('fusion_gene_aligned_pairs.rds')
se <- readRDS('rnaseq_read_count_entrez_id.rds')
molm13.fpm <- 1e6 * elementNROWS(l) / colSums(assay(se))

condition <- c(
    rep.int( c('ADA', 'DCD'), times =  rep(length(molm13.fpm)/2,2) ),
    rep.int( c('ADA', 'DCD'), times =  rep(length(lsk.fpm)/2,2) ),
    rep.int( c('ADA', 'DCD'), times =  rep(length(hsc.fpm)/2,2) ),
    rep.int( c('ADA', 'DCD'), times =  rep(length(lsc.fpm)/2,2) ) )
condition <- factor( condition, levels =  c('ADA', 'DCD') )
cell <- rep.int( c( "MOLM13", "LSK", "HSC", "LSC" ), 
                times = c(length(molm13.fpm), length(lsk.fpm), length(hsc.fpm), length(lsc.fpm) ) )
cell <- factor( cell, levels = c( "MOLM13", "LSK", "HSC", "LSC" ) )
df <- data.frame( fpm = c( molm13.fpm, lsk.fpm, hsc.fpm, lsc.fpm ),
           cell = cell, condition = condition )
ggplot( df, aes( x = cell, y = fpm, color = condition, order = rev( order( condition ) )  ) ) + 
geom_point( position = position_jitter(width = .1) ) + scale_color_manual(values=c('grey25','grey75')) +
xlab("Cell") + ylab("FPM") + theme_classic()
