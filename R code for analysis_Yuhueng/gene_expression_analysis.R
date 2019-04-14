library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DESeq2)
library(stringr)
library(BiocParallel)
library(openxlsx)

se <- readRDS("rnaseq_read_count_entrez_id.rds")

df.list <- list()
compartments <- c( "MPP2", "MPP4", "LT", "ST" )
for( i in compartments ){
se1 <- se[,colData(se)$cell == i]
fpkm.mat <- fpkm(DESeqDataSet(se1,~condition))
fpkm.mat <- data.frame(
ADA.fpkm = rowMeans(fpkm.mat[,colData(se1)$condition == 'ADA']),
DCD.fpkm = rowMeans(fpkm.mat[,colData(se1)$condition == 'DCD']),
MIG.fpkm = rowMeans(data.frame(fpkm.mat[,colData(se1)$condition == 'MIG']) )
                      )
rowmask <- fpkm.mat$ADA.fpkm >= 1 |
fpkm.mat$DCD.fpkm >= 1 |
fpkm.mat$MIG.fpkm >= 1

df.list[[i]] <- data.frame( 
                           entrez.id = rownames(se1)[rowmask], 
                           gene.symbol = unlist( mget( rownames(se1)[rowmask], org.Mm.egSYMBOL, ifnotfound = NA ) ), 
                           fpkm.mat[rowmask,], assay(se1)[rowmask,] )
}

write.xlsx( df.list, file = "mouse_hsc_mettl3_gene_expression.xlsx", row.names = F )
