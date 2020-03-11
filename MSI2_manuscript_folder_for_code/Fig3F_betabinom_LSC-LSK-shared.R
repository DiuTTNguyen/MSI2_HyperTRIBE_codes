library(openxlsx)
library(parallel)
library(stringr)
library(bbmle)
library(VGAM)
library(DEXSeq)
library(dplyr)

# setwd("/Users/karen/mount/chuk/Fig3B_betabinomial_on_LSC_and_LSK_unique_targets/")

df <- read.csv( "lsc.lsk.combined.betabinom.input_diff.freq_greaterthanequal_0.6.csv" )

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

adar.index <- which( grepl("ADA", colnames(df)) )
df1 <- df[,adar.index]
df1 <- cbind(df[,1:9], df1)

ref.counts <- df1[,grep('ref.count', colnames(df1))]
alt.counts <- df1[,grep('alt.count', colnames(df1))]

lsc.adar.index <- which( grepl("ADA.*DsRed", colnames(ref.counts)) )
lsk.adar.index <- which( grepl("ADA.A_IGO|ADA.B_IGO|ADA.C_IGO", colnames(ref.counts)) )

# lsc.dcd.index <- which( grepl("DCD.*DsRed|DCD-B-Ds-Red", colnames(ref.counts)) )
# lsk.dcd.index <- which( grepl("DCD-A_IGO|DCD-B_IGO|DCD-C_IGO", colnames(ref.counts)) )
# 
# lsc.mig.index <- which( grepl("MIG.*DsRed", colnames(ref.counts)) )
# lsk.mig.index <- which( grepl("MIG-A_IGO|MIG-B_IGO|MIG-C_IGO", colnames(ref.counts)) )

lrt.res <- mclapply( 1:nrow(df1), function( i ){ 
       counts <- as.integer( df1[i,-(1:9)] )
       ref <- counts[ seq( 1, length( counts ), 2 ) ]
       alt <- counts[ seq( 2, length( counts ), 2 ) ]
       fit0 <- mle.custom.h0( ref[ c(lsc.adar.index, lsk.adar.index) ], alt[ c(lsc.adar.index, lsk.adar.index) ] )
       fit1 <- mle.custom.h1( ref[ lsc.adar.index ], alt[ lsc.adar.index ], ref[ lsk.adar.index ], alt[ lsk.adar.index ] )
       p.value <- pchisq( 2 *( fit0@min - fit1@min ), 1, lower.tail=F )
       value <- fit1@min
       code <- fit1@details$convergence
       return( c( p.value = p.value, value = value, code = code ) )
}, mc.cores = 12 )
lrt.res <- data.frame( do.call( 'rbind', lrt.res ) )
saveRDS(lrt.res, "lrt.res_LSC_LSK_combined_diff.freq_greaterthanequal_0.6.RDS")

# Additional columns
ref.counts <- df[,grep('ref.count', colnames(df))]
alt.counts <- df[,grep('alt.count', colnames(df))]

lsc.adar.index <- which( grepl("ADA.*DsRed", colnames(ref.counts)) )
lsk.adar.index <- which( grepl("ADA.A_IGO|ADA.B_IGO|ADA.C_IGO", colnames(ref.counts)) )

lsc.dcd.index <- which( grepl("DCD.*DsRed|DCD.B.Ds.Red", colnames(ref.counts)) )
lsk.dcd.index <- which( grepl("DCD.A_IGO|DCD.B_IGO|DCD.C_IGO", colnames(ref.counts)) )
 
lsc.mig.index <- which( grepl("MIG.*DsRed", colnames(ref.counts)) )
lsk.mig.index <- which( grepl("MIG.A_IGO|MIG.B_IGO|MIG.C_IGO", colnames(ref.counts)) )

alt.freq <- alt.counts / ( ref.counts + alt.counts )
df.stats <- data.frame(diff.frequency = rowMeans(alt.freq[, lsc.adar.index ], na.rm=T) - rowMeans(alt.freq[,lsk.adar.index],na.rm=T),
                       LSC.diff.frequency = rowMeans(alt.freq[, lsc.adar.index ], na.rm=T) - rowMeans(alt.freq[,c(lsc.dcd.index, lsc.mig.index)],na.rm=T),
                       LSC.ADA.frequency = rowMeans(alt.freq[,lsc.adar.index], na.rm=T),
                       LSC.DCD.frequency = rowMeans(alt.freq[,lsc.dcd.index], na.rm=T),
                       LSC.MIG.frequency = rowMeans(alt.freq[,lsc.mig.index], na.rm=T),
                       
                       LSK.diff.frequency = rowMeans(alt.freq[, lsk.adar.index ], na.rm=T) - rowMeans(alt.freq[,c(lsk.dcd.index, lsk.mig.index)],na.rm=T),
                       LSK.ADA.frequency = rowMeans(alt.freq[,lsk.adar.index], na.rm=T),
                       LSK.DCD.frequency = rowMeans(alt.freq[,lsk.dcd.index], na.rm=T),
                       LSK.MIG.frequency = rowMeans(alt.freq[,lsk.mig.index], na.rm=T) )
df.stats <- cbind( df.stats, lrt.res )
write.csv(df.stats, "df.stats_LSC_LSK_combined_diff.freq_greaterthanequal_0.6.csv", row.names = FALSE)

rowmask <- rowSums( alt.counts[c(lsc.adar.index, lsk.adar.index)] != 0 ) > 1 & # TRUE/FALSE on which row has sum of alt.counts > 1 (Karen)
( rowMeans( ref.counts[ lsc.adar.index ] ) + rowMeans( alt.counts[ lsc.adar.index ] ) >= 5 ) & # Seems to be ADAR samples for this line (Karen)
( rowMeans( ref.counts[ lsk.adar.index ] ) + rowMeans( alt.counts[ lsk.adar.index ] ) >= 5 ) &
( rowMeans( ref.counts[c(lsc.adar.index, lsk.adar.index)] ) + rowMeans( alt.counts[c(lsc.adar.index, lsk.adar.index)] ) >= 5 ) &
lrt.res$value > 0 & 
!is.na(df$gene.symbol) 

df.stats <- df.stats[ rowmask, ]
df.stats$p.adj <- p.adjust( df.stats$p.value, 'BH' )
df.stats$value <- NULL
df.stats$code <- NULL
df1 <- cbind( df[ rowmask, c(1:9) ], df.stats, df[ rowmask, -(1:9) ] )
write.csv(df1, "Mouse_LSC_LSK_ADA_DCD_MIG_snps_count_dedupped_KNOWN_SNPs_removed_LSC_LSK_combined_diff.freq_greaterthanequal_0.6_significance.csv", row.names = FALSE)

# Add gene FPKM and # edits
gene.num.edits <- table( df1$entrez.id )
se <- readRDS( "Yuheng_files/rnaseq_read_count_entrez_id.rds" )
fpkm.mat <- fpkm( DESeqDataSet( se, ~1 ) ) # originally ~condition but R gives error cause it was never defined -Karen

lsc.adar.index <- which( grepl("ADA.*DsRed", colnames(fpkm.mat)) )
lsk.adar.index <- which( grepl("ADA-A_IGO|ADA-B_IGO|ADA-C_IGO", colnames(fpkm.mat)) )

lsc.dcd.index <- which( grepl("DCD.*DsRed|DCD-B-Ds-Red", colnames(fpkm.mat)) )
lsk.dcd.index <- which( grepl("DCD-A_IGO|DCD-B_IGO|DCD-C_IGO", colnames(fpkm.mat)) )

lsc.mig.index <- which( grepl("MIG.*DsRed", colnames(fpkm.mat)) )
lsk.mig.index <- which( grepl("MIG-A_IGO|MIG-B_IGO|MIG-C_IGO", colnames(fpkm.mat)) )

fpkm.mat <- data.frame( LSC.ADAR.fpkm = rowMeans( fpkm.mat[,lsc.adar.index]),
                        LSC.DCD.fpkm = rowMeans( fpkm.mat[,lsc.dcd.index]),
                        LSC.MIG.fpkm = rowMeans( fpkm.mat[,lsc.mig.index]),
                        LSK.ADAR.fpkm = rowMeans( fpkm.mat[,lsk.adar.index]),
                        LSK.DCD.fpkm = rowMeans( fpkm.mat[,lsk.dcd.index]),
                        LSK.MIG.fpkm = rowMeans( fpkm.mat[,lsk.mig.index]) )
fpkm.mat.merge <- data.frame( entrez.id = names(gene.num.edits), gene.num.edits = as.integer(gene.num.edits), fpkm.mat[ names( gene.num.edits ), ] )
df1 <- merge( df1, fpkm.mat.merge, by = 'entrez.id' )
write.csv(df1, "Mouse_LSC_LSK_ADA_DCD_MIG_snps_count_dedupped_KNOWN_SNPs_removed_LSC_LSK_combined_diff.freq_greaterthanequal_0.6_significance_fpkm.csv", row.names = FALSE)


