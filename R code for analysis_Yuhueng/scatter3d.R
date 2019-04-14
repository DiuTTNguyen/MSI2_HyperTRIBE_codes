library(openxlsx)
library(scatterplot3d)

df <- read.xlsx( "mouse_lsk_snp_counts_dedupped.xlsx" )

rowmask <- ( df$`Sample_Hyper-dADAR_12617_IGO_08269_2.ref.count` + df$`Sample_Hyper-dADAR_12617_IGO_08269_2.alt.count` >= 25 ) &
( df$`Sample_Hyper-dADAR_121417_IGO_08269_5.ref.count` + df$`Sample_Hyper-dADAR_121417_IGO_08269_5.alt.count` >= 25 ) &
( df$`Sample_Hyper-dADAR_121817_IGO_08269_8.ref.count` + df$`Sample_Hyper-dADAR_121817_IGO_08269_8.alt.count` >= 25 )
num.sites <- sum( rowmask )
filler <- rep( 0, num.sites )

x <- df$`Sample_Hyper-dADAR_12617_IGO_08269_2.alt.count` / ( df$`Sample_Hyper-dADAR_12617_IGO_08269_2.ref.count` + df$`Sample_Hyper-dADAR_12617_IGO_08269_2.alt.count` )
y <- df$`Sample_Hyper-dADAR_121417_IGO_08269_5.alt.count` / ( df$`Sample_Hyper-dADAR_121417_IGO_08269_5.ref.count` + df$`Sample_Hyper-dADAR_121417_IGO_08269_5.alt.count` )
z <- df$`Sample_Hyper-dADAR_121817_IGO_08269_8.alt.count` / ( df$`Sample_Hyper-dADAR_121817_IGO_08269_8.ref.count` + df$`Sample_Hyper-dADAR_121817_IGO_08269_8.alt.count` )

cor.xy <- cor( x[ rowmask ], y[ rowmask ] )
cor.xz <- cor( x[ rowmask ], z[ rowmask ] )
cor.yz <- cor( y[ rowmask ], z[ rowmask ] )
sp3d <- scatterplot3d( x[ rowmask ], z[ rowmask ],  y[ rowmask ],  
                      xlab = "dADAR-1", ylab = "dADAR-3", zlab = "dADAR-2",
              grid = F,  angle = 45, pch = 3, color = 'lightgrey' )
sp3d$points3d( filler, z[ rowmask ], y[ rowmask ], pch = 19, col = 'grey60')
sp3d$points3d( x[ rowmask ], z[ rowmask ], filler, pch = 19, col = 'grey45' )
sp3d$points3d( x[ rowmask ], filler, y[ rowmask ], pch = 19, col = 'grey30')
text( sp3d$xyz.convert( .8, 1, .075 ), labels = substitute( paste( r[13], " = ", cor ), list( cor = sprintf( "%.3f", cor.xz ) ) ) )
text( sp3d$xyz.convert( .55, 0, .875 ), labels = substitute( paste( r[12], " = ", cor ), list( cor = sprintf( "%.3f", cor.xy ) ) ) )
text( sp3d$xyz.convert( .15, 1, .8 ), labels = substitute( paste( r[23], " = ", cor ), list( cor = sprintf( "%.3f", cor.yz ) ) ) )
