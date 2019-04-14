library(openxlsx)
library(parallel)
library(stringr)

library(bbmle)
library(VGAM)

library(DEXSeq)

df <- read.xlsx( "mouse_lsk_snp_counts_dedupped.xlsx", sheet = 1 )

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

lrt.res <- mclapply( 1:nrow(df), function( i ){ 
       counts <- as.integer( df[i,-(1:8)] )
       ref <- counts[ seq( 1, length( counts ), 2 ) ]
       alt <- counts[ seq( 2, length( counts ), 2 ) ]
       fit0 <- mle.custom.h0( ref, alt )
       fit1 <- mle.custom.h1( ref[1:2], alt[1:2], ref[-(1:2)], alt[-(1:2)] )
       p.value <- pchisq( 2 *( fit0@min - fit1@min ), 1, lower.tail=F )
       value <- fit1@min
       code <- fit1@details$convergence
       return( c( p.value = p.value, value = value, code = code ) )
}, mc.cores = 12 )
lrt.res <- data.frame( do.call( 'rbind', lrt.res ) )

# Additional columns
# p.value <- unlist( tmp )
alt.freq <- alt.counts / ( ref.counts + alt.counts )
df.stats <- data.frame( 
                       diff.frequency = rowMeans(alt.freq[,c(1,2,3)],na.rm=T) - rowMeans(alt.freq[,-c(1,2,3)],na.rm=T),
                       ADA.frequency = rowMeans( alt.freq[,c(1,2,3)],na.rm=T ),
                       DCD.frequency = rowMeans( alt.freq[,c(4,5,6)],na.rm=T ),
                       MIG.frequency = rowMeans( alt.freq[,c(7,8,9)],na.rm=T ) )

df.stats <- cbind( df.stats, lrt.res )

rowmask <- rowSums( alt.counts != 0 ) > 1 &
( rowMeans( ref.counts[1:2] ) + rowMeans( alt.counts[1:2] ) >= 5 ) &
( rowMeans( ref.counts ) + rowMeans( alt.counts ) >= 5 ) &
lrt.res$value > 0 & 
!is.na(df$gene.symbol) 

df.stats <- df.stats[ rowmask, ]
df.stats$p.adj <- p.adjust( df.stats$p.value, 'BH' )
df.stats$value <- NULL
df.stats$code <- NULL
df1 <- cbind( df[ rowmask, 1:8 ], df.stats, df[ rowmask, -(1:8) ] )
df1 <- df1[ which( df1$diff.frequency > 0 & df1$p.adj < .05 ), ]
                       


# Transform
# c( mu = fit@coef[1] / ( fit@coef[1] + fit@coef[2] ), 
# rho = 1 / ( fit@coef[1] + fit@coef[2] + 1 ) )
# VGAM implementation
# Coef( vglm( cbind( alt, ref ) ~ 1, betabinomial ) )

# Issues: 
# 1. rho ~ 0 in most cases; Binomial approximation?
# 2. Alt read all-zero in many cases; Pseudo count necessary?

rowmask <- apply( ref.counts, 1, function(x) !all( x== 0))
df1 <- df[ rowmask, ]
# Comparison with DEXSeq
# Reshape the count matrix
idx <- as.vector( rbind( 1:nrow(df1), (1:nrow(df1))+nrow(df1)) )
all.counts <- rbind( as.matrix(ref.counts), as.matrix(alt.counts) )[ idx, ] 
colnames( all.counts ) <- str_split_fixed( colnames( all.counts ), '\\.', 2 )[,1]
condition <- str_split_fixed( colnames( ref.counts ), "_", 3 )[,2]
dexds <- DEXSeqDataSet( countData = all.counts, 
                       sampleData = data.frame( condition = condition ), 
                       featureID = rep( c('ref', 'alt'), nrow( df1 ) ),
                       groupID = rep( paste( df1[,1], df1[,2], sep = '_' ), 2 ), 
                       design =  ~ sample + exon + condition:exon )

dexds <- estimateSizeFactors(dexds)
dexds <- estimateDispersions(dexds, BPPARAM = MulticoreParam( workers = 20 ))
dexds <- testForDEU( dexds )
dexds <- estimateExonFoldChanges( dexds, fitExpToVar="condition", BPPARAM = MulticoreParam( workers = 20 ) )
dex.res <-  DEXSeqResults( dexds )

# Filtering

# 1. Only 1 non-zero alt / Low read counts
# 2. Calculate FDR
# 3. FDR < 0.01 & diff > 0 & Convergence
df.stats <- cbind( df.stats, lrt.res )

rowmask <- rowSums( alt.counts != 0 ) > 1 &
( rowSums( ref.counts[1:3] ) + rowSums( alt.counts[1:3] ) >= 15 ) &
( rowSums( ref.counts ) + rowSums( alt.counts ) >= 45 ) &
lrt.res$value > 0 & 
!is.na(df$gene.symbol) 

df.stats <- df.stats[ rowmask, ]
df.stats$p.adj <- p.adjust( df.stats$p.value, 'BH' )
df.stats$value <- NULL
df.stats$code <- NULL
df1 <- cbind( df[ rowmask, 1:8 ], df.stats, df[ rowmask, -(1:8) ] )
df1 <- df1[ which( df1$diff.frequency > 0 & df1$p.adj < .01 ), ]

# Add gene FPKM and # edits
gene.num.edits <- table( df1$entrez.id )
se <- readRDS( "rnaseq_read_count_entrez_id.rds" )
fpkm.mat <- fpkm( DESeqDataSet( se, ~condition ) )
fpkm.mat <- data.frame( ADA.fpkm = rowMeans( fpkm.mat[,1:3]),
                       DCD.fpkm = rowMeans( fpkm.mat[,4:6]),
                       MIG.fpkm = rowMeans( fpkm.mat[,7:9]) )
fpkm.mat <- data.frame( entrez.id = names(gene.num.edits), gene.num.edits = as.integer(gene.num.edits), fpkm.mat[ names( gene.num.edits ), ] )
df1 <- merge( df1, fpkm.mat, by = 'entrez.id' )
