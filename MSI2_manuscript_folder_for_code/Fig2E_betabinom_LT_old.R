library(openxlsx)
library(parallel)
library(stringr)
library(bbmle)
library(VGAM)
library(DEXSeq)
library(dplyr)

setwd("/data/leslie/chuk/LT_vs_MPPs/")

df1 <- read.csv( "prepare_betabinom_input/Mouse_HSPC_LT-unique_betabinom_input.csv" )
df1 <- df1 %>% select(-"X")

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

lt.adar.index <- which( grepl("ADA.*LT", colnames(ref.counts)) )
# lt.mig.index <- which( grepl("MIG.*LT", colnames(ref.counts)) )
# lt.dcd.index <- which( grepl("DCD.*LT", colnames(ref.counts)) )

st.adar.index <- which( grepl("ADA.*ST", colnames(ref.counts)) )
# st.mig.index <- which( grepl("MIG.*ST", colnames(ref.counts)) )
# st.dcd.index <- which( grepl("DCD.*ST", colnames(ref.counts)) )

mpp2.adar.index <- which( grepl("ADA.*MPP2", colnames(ref.counts)) )
# mpp2.mig.index <- which( grepl("MIG.*MPP2", colnames(ref.counts)) )
# mpp2.dcd.index <- which( grepl("DCD.*MPP2", colnames(ref.counts)) )

mpp4.adar.index <- which( grepl("ADA.*MPP4", colnames(ref.counts)) )
# mpp4.mig.index <- which( grepl("MIG.*MPP4", colnames(ref.counts)) )
# mpp4.dcd.index <- which( grepl("DCD.*MPP4", colnames(ref.counts)) )

lrt.res <- mclapply( 1:nrow(df1), function( i ){
       counts <- as.integer( df1[i,-(1:9)] )
       ref <- counts[ seq( 1, length( counts ), 2 ) ]
       alt <- counts[ seq( 2, length( counts ), 2 ) ]
       fit0 <- mle.custom.h0( ref, alt )
       fit1 <- mle.custom.h1( ref[ lt.adar.index ], alt[ lt.adar.index ],
                              ref[ -lt.adar.index ], alt[ -lt.adar.index ] )
       p.value <- pchisq( 2 *( fit0@min - fit1@min ), 1, lower.tail=F )
       value <- fit1@min
       code <- fit1@details$convergence
       return( c( p.value = p.value, value = value, code = code ) )
}, mc.cores = 12 )
lrt.res <- data.frame( do.call( 'rbind', lrt.res ) )
saveRDS(lrt.res, "betabinom/lrt.res_HSPC_LT-unique_vs_others.RDS")

# Additional columns
alt.freq <- alt.counts / ( ref.counts + alt.counts )
df.stats <- data.frame(diff.frequency = rowMeans(alt.freq[, lt.adar.index ], na.rm=T) -
                         rowMeans(alt.freq[,-lt.adar.index],na.rm=T),

                       # LT.diff.frequency = rowMeans(alt.freq[, lt.adar.index ], na.rm=T) -
                       #   rowMeans(alt.freq[,c(lt.mig.index, lt.dcd.index)],na.rm=T),
                       LT.ADA.frequency = rowMeans(alt.freq[, lt.adar.index ], na.rm=T),
                       # LT.MIG.frequency = rowMeans(alt.freq[, lt.mig.index ], na.rm=T),
                       # LT.DCD.frequency = rowMeans(alt.freq[, lt.dcd.index ], na.rm=T),

                       # ST.diff.frequency = rowMeans(alt.freq[, st.adar.index ], na.rm=T) -
                       #   rowMeans(alt.freq[,c(st.mig.index, st.dcd.index)],na.rm=T),
                       ST.ADA.frequency = rowMeans(alt.freq[, st.adar.index ], na.rm=T),
                       # ST.MIG.frequency = alt.freq[, st.mig.index ],
                       # ST.DCD.frequency = rowMeans(alt.freq[, st.dcd.index ], na.rm=T),

                       # MPP2.diff.frequency = rowMeans(alt.freq[, mpp2.adar.index ], na.rm=T) -
                       #   rowMeans(alt.freq[,c(mpp2.mig.index, mpp2.dcd.index)],na.rm=T),
                       MPP2.ADA.frequency = rowMeans(alt.freq[, mpp2.adar.index ], na.rm=T),
                       # MPP2.MIG.frequency = rowMeans(alt.freq[, mpp2.mig.index ], na.rm=T),
                       # MPP2.DCD.frequency = rowMeans(alt.freq[, mpp2.dcd.index ], na.rm=T),

                       # MPP4.diff.frequency = rowMeans(alt.freq[, mpp4.adar.index ], na.rm=T) -
                       #   rowMeans(alt.freq[,c(mpp4.mig.index, mpp4.dcd.index)],na.rm=T),
                       MPP4.ADA.frequency = rowMeans(alt.freq[, mpp4.adar.index ], na.rm=T)
                       # MPP4.MIG.frequency = rowMeans(alt.freq[, mpp4.mig.index ], na.rm=T),
                       # MPP4.DCD.frequency = rowMeans(alt.freq[, mpp4.dcd.index ], na.rm=T) 
                       )
df.stats <- cbind( df.stats, lrt.res )
write.csv(df.stats, "betabinom/df.stats_HSPC_LT-unique_vs_others.csv", row.names = FALSE)

rowmask <- rowSums( alt.counts != 0 ) > 1 & # TRUE/FALSE on which row has sum of alt.counts > 1 (Karen)
( rowMeans( ref.counts[ lt.adar.index ] ) + rowMeans( alt.counts[ lt.adar.index ] ) >= 5 ) & 
( rowMeans( ref.counts ) + rowMeans( alt.counts ) >= 5 ) &
lrt.res$value > 0 &
!is.na(df1$gene.symbol)

df.stats <- df.stats[ rowmask, ]
df.stats$p.adj <- p.adjust( df.stats$p.value, 'BH' )
df.stats$value <- NULL
df.stats$code <- NULL
df1 <- cbind( df1[ rowmask, c(1:9) ], df.stats, df1[ rowmask, -(1:9) ] )
write.csv(df1, "betabinom/Mouse_HSPC_LT-unique_vs_others_snp_counts_significance.csv", row.names = FALSE)

# Add gene FPKM and # edits
df1 <- read.csv("betabinom/Mouse_HSPC_LT-unique_vs_others_snp_counts_significance.csv")
gene.num.edits <- table( df1$entrez.id )
se <- readRDS( "snp_counts/rnaseq_read_count_entrez_id.rds" )
fpkm.mat <- fpkm( DESeqDataSet( se, ~1 ) ) # originally ~condition but R gives error cause it was never defined -Karen

lt.adar.index <- which( grepl("ADA.*LT", colnames(fpkm.mat)) )
# lt.mig.index <- which( grepl("MIG.*LT", colnames(fpkm.mat)) )
# lt.dcd.index <- which( grepl("DCD.*LT", colnames(fpkm.mat)) )

st.adar.index <- which( grepl("ADA.*ST", colnames(fpkm.mat)) )
# st.mig.index <- which( grepl("MIG.*ST", colnames(fpkm.mat)) )
# st.dcd.index <- which( grepl("DCD.*ST", colnames(fpkm.mat)) )

mpp2.adar.index <- which( grepl("ADA.*MPP2", colnames(fpkm.mat)) )
# mpp2.mig.index <- which( grepl("MIG.*MPP2", colnames(fpkm.mat)) )
# mpp2.dcd.index <- which( grepl("DCD.*MPP2", colnames(fpkm.mat)) )

mpp4.adar.index <- which( grepl("ADA.*MPP4", colnames(fpkm.mat)) )
# mpp4.mig.index <- which( grepl("MIG.*MPP4", colnames(fpkm.mat)) )
# mpp4.dcd.index <- which( grepl("DCD.*MPP4", colnames(fpkm.mat)) )


fpkm.mat <- data.frame( LT.ADA.fpkm = rowMeans( fpkm.mat[,lt.adar.index]),
                        
                        ST.ADA.fpkm = rowMeans( fpkm.mat[,st.adar.index]),
                        
                        MPP2.ADA.fpkm = rowMeans( fpkm.mat[,mpp2.adar.index]),
                        
                        MPP4.ADA.fpkm = rowMeans( fpkm.mat[,mpp4.adar.index]) )

fpkm.mat.merge <- data.frame( entrez.id = names(gene.num.edits), gene.num.edits = as.integer(gene.num.edits), fpkm.mat[ names( gene.num.edits ), ] )
df1 <- merge( df1, fpkm.mat.merge, by = 'entrez.id' )
write.csv(df1, "betabinom/Mouse_HSPC_LT-unique_vs_others_snp_counts_significance_fpkm.csv", row.names = FALSE)


