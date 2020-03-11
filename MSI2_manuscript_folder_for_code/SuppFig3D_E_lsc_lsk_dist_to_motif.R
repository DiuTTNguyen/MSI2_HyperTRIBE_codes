# LSC and LSK distance to motif probability density function plots
# Karen Chu
# 2/12/2019

# Purpose: Plotting probability density function of Msi2 HyperTRIBE edited sites to nearest MSI2 motif in LSC and LSK.

library(ggplot2)

setwd("/Users/chuk/Documents/MSI2-hyperTRIBE/LSC_LSK_distance_to_nearest_motif/")

lsc <- read.csv("data/lsc_with_distance.csv")
lsk <- read.csv("data/lsk_with_distance.csv")

# Filter by fpkm >= 5 and diff.frequency >= 0.01
lsc <- lsc [ lsc$ADA.fpkm >= 5 & lsc$DCD.fpkm >= 5 & lsc$MIG.fpkm >=5 & lsc$diff.frequency >= 0.1, ]
lsk <- lsk [ lsk$ADA.fpkm >= 5 & lsk$DCD.fpkm >= 5 & lsk$MIG.fpkm >=5 & lsk$diff.frequency >= 0.1, ]

# Plot LSC PDF plot: lsc_dist_to_motif.png
pdf("figures/lsc_dist_to_motif.pdf", 15, 12, useDingbats = F )
ggplot(lsc, aes(x=lsc$dist.to.motif)) + geom_density(fill="dodgerblue") +
  scale_x_continuous(limits = c(-500, 500)) +
  labs(title="LSC Distance to nearest motif\n") +
  labs(x="\nDistance (bp)") +
  labs(y="Density\n") +
  theme_classic() +
  theme(plot.title = element_text(size=40), 
        axis.text=element_text(size=40), 
        axis.text.x = element_text(angle = 90),
        axis.title.x = element_text(size=40), 
        axis.title.y = element_text(size=40))
dev.off()

# Plot LSK PDF plot: lsk_dist_to_motif.png
pdf("figures/lsk_dist_to_motif.pdf", 15, 12, useDingbats = F )
ggplot(lsk, aes(x=lsk$dist.to.motif)) + geom_density(fill="dodgerblue") +
  scale_x_continuous(limits = c(-500, 500)) +
  labs(title="LSK Distance to nearest motif\n") +
  labs(x="\nDistance (bp)") +
  labs(y="Density\n") +
  theme_classic() +
  theme(plot.title = element_text(size=40), 
        axis.text=element_text(size=40), 
        axis.text.x = element_text(angle = 90),
        axis.title.x = element_text(size=40), 
        axis.title.y = element_text(size=40))
dev.off()

# Histogram
ggplot(data=lsc, aes(lsc$dist.to.motif)) + 
  geom_histogram(breaks=seq(-500,500, by=50), 
                 col="white", 
                 fill="blue", 
                 alpha=0.2) + 
  labs(title="LSC Distance to nearest motif") +
  labs(x="Distance (bp)") +
  labs(y="Count")
