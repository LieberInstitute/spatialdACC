setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("SpatialExperiment")
library("patchwork")
library("tidyverse")
library("viridis")
library("pheatmap")
library("ComplexHeatmap")
library("scater")
library("bluster")
library("sessioninfo")
library("here")
library("schex")
library("svglite")
library("dplyr")
library("patchwork")


load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata"))
df_2 <- as.data.frame(colData(spe))
df_2 <- df_2 %>%
    group_by(brnum) %>%
    count(layer)

p2 <- ggplot(data = df_2, aes(x=layer, y=n, fill=brnum)) +
    geom_bar(stat="identity") +
    ylab("Number of Spots") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    xlab("Spatial Domain") +
    theme(legend.position="none")


png(file = here::here("plots", "08_clustering", "spatial_domain_barplot.png"), height = 6, width = 6, unit="in",res=300)
p2
dev.off()
