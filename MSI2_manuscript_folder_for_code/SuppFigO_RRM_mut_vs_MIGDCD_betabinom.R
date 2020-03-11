library(openxlsx)
library(parallel)
library(stringr)
library(bbmle)
library(VGAM)
library(DEXSeq)
library(dplyr)

# setwd("/Users/karen/mount/chuk/MSI2_RRM_deletion_and_5_mutations/Project_08269_F/5mut_vs_MIG/")

df1 <- read.xlsx( "MOLM13_MSI2_5_mutations_vs_MIG_snp_counts.xlsx" )
df1 <- subset(df1, dbsnp=="")

# MLE 

mle.custom.h0 <- function( ref, alt, debug = F ){
    x <- alt
    size <- ref + alt

    bbll <- function( prob, rho ){
        if( ( prob > 0 ) & ( rho > 0 ) & ( prob < 1 ) & ( rho < 1 ) )
            -sum( dbetabinom( x, size, prob, rho, log = T ) )
        else NA 
    }

    fit <- mle2( bbll, # nobs = length(x), 
               # start = list( shape1 = 10, shape2 = 5 ),
               start = list( prob = .1, rho = .5 ),
               method = "Nelder-Mead", # hessian = FALSE, 
               control = list( maxit = 1e5, trace = as.integer(debug) ) ) 
    return( fit )
}

mle.custom.h1 <- function( ref1, alt1, ref2, alt2, debug = F ){
    x1 <- alt1
    size1 <- ref1 + alt1
    prob1.init <- mean(x1/size1)  + 1e-3
    if( is.na(prob1.init) | prob1.init >= 1 | prob1.init <= 0 ){
        prob1.init <- .05
    } 
    x2 <- alt2
    size2 <- ref2 + alt2
    prob2.init <- mean(x2/size2)  + 1e-3

    if( is.na(prob2.init) | prob2.init >= 1 | prob2.init <= 0 ){
        prob2.init <- .05
    } 

    bbll <- function( prob1, prob2, rho ){
        if( ( prob1 > 0 ) & ( prob2 > 0 ) & ( rho > 0 ) &
            ( prob1 < 1 ) & ( prob2 < 1 ) & ( rho < 1 ) )
            -( sum( dbetabinom( x1, size1, prob1, rho, log = T ) ) + 
               sum( dbetabinom( x2, size2, prob2, rho, log = T ) ) )
        else NA 
    }

    fit <- mle2( bbll,
                start = list( prob1 = prob1.init, prob2 = prob2.init, rho = .1 ),
               method = "Nelder-Mead", # hessian = FALSE, 
               control = list( maxit = 1e5, trace = as.integer(debug) ) ) 
    return( fit )
}

ref.counts <- df1[,grep('ref.count', colnames(df1))]
alt.counts <- df1[,grep('alt.count', colnames(df1))]

mut.index <- which( grepl("5mut", colnames(ref.counts)) )
dcd.index <- which( grepl("DCD", colnames(ref.counts)) )
mig.index <- which( grepl("MIG", colnames(ref.counts)) )

lrt.res <- mclapply( 1:nrow(df1), function( i ){ 
       counts <- as.integer( df1[i,-(1:9)] )
       ref <- counts[ seq( 1, length( counts ), 2 ) ]
       alt <- counts[ seq( 2, length( counts ), 2 ) ]
       fit0 <- mle.custom.h0( ref, alt )
       fit1 <- mle.custom.h1( ref[ mut.index ], alt[ mut.index ], 
                              ref[ c(dcd.index, mig.index) ], alt[ c(dcd.index, mig.index) ] )
       p.value <- pchisq( 2 *( fit0@min - fit1@min ), 1, lower.tail=F )
       value <- fit1@min
       code <- fit1@details$convergence
       return( c( p.value = p.value, value = value, code = code ) )
}, mc.cores = 12 )
lrt.res <- data.frame( do.call( 'rbind', lrt.res ) )
saveRDS(lrt.res, "lrt.res.RDS")

# Additional columns
alt.freq <- alt.counts / ( ref.counts + alt.counts )
df.stats <- data.frame(diff.frequency = rowMeans(alt.freq[, mut.index ], na.rm=T) - 
                         rowMeans(alt.freq[,c(dcd.index,mig.index)],na.rm=T),
                       mut.frequency = rowMeans(alt.freq[,mut.index], na.rm=T),
                       dcd.frequency = rowMeans(alt.freq[,dcd.index], na.rm=T),
                       mig.frequency = rowMeans(alt.freq[,mig.index], na.rm=T) )
df.stats <- cbind( df.stats, lrt.res )
write.csv(df.stats, "df.stats.csv", row.names = FALSE)

rowmask <- rowSums( alt.counts != 0 ) > 1 & # TRUE/FALSE on which row has sum of alt.counts > 1 (Karen)
( rowMeans( ref.counts[ mut.index ] ) + rowMeans( alt.counts[ mut.index ] ) >= 5 ) & # Seems to be ADAR samples for this line (Karen)
( rowMeans( ref.counts ) + rowMeans( alt.counts ) >= 5 ) &
lrt.res$value > 0 & 
!is.na(df1$gene.symbol) 

df.stats <- df.stats[ rowmask, ]
df.stats$p.adj <- p.adjust( df.stats$p.value, 'BH' )
df.stats$value <- NULL
df.stats$code <- NULL
df1 <- cbind( df1[ rowmask, c(1:9) ], df.stats, df1[ rowmask, -(1:9) ] )
write.csv(df1, "MOLM13_5mut_vs_MIGDCD_snp_counts_significance.csv", row.names = FALSE)

# Add gene FPKM and # edits
gene.num.edits <- table( df1$entrez.id )
se <- readRDS( "rnaseq_read_count_entrez_id.rds" )
fpkm.mat <- fpkm( DESeqDataSet( se, ~1 ) ) # originally ~condition but R gives error cause it was never defined -Karen

mut.index <- which( grepl("5mut", colnames(fpkm.mat)) )
dcd.index <- which( grepl("DCD", colnames(fpkm.mat)) )
mig.index <- which( grepl("MIG", colnames(fpkm.mat)) )

fpkm.mat <- data.frame( mut.fpkm = rowMeans( fpkm.mat[,mut.index]),
                        dcd.fpkm = rowMeans( fpkm.mat[,dcd.index]),
                        mig.fpkm = rowMeans( fpkm.mat[,mig.index]) )
fpkm.mat.merge <- data.frame( entrez.id = names(gene.num.edits), gene.num.edits = as.integer(gene.num.edits), fpkm.mat[ names( gene.num.edits ), ] )
df1 <- merge( df1, fpkm.mat.merge, by = 'entrez.id' )
write.csv(df1, "MOLM13_5mut_vs_MIGDCD_snp_counts_significance_fpkm.csv", row.names = FALSE)


