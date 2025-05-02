setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggspavis)
library(gridExtra)
library(ggrepel)
library(here)
library(scran)
library(scry)

load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9.Rdata"))

# compare library size of WM1, WM2, WM-CC
# overall and separately by sample
# make a dot plot of the library size of each sample

# create a dataframe by taking the average of colData(spe)$sizeFactor for each sample in colData(spe)$sample_id and cluster in colData(spe)$PRECAST_cluster

# Filter for relevant clusters
wm_clusters <- c("3", "8", "7")

# Create a dataframe with the average library size (sizeFactor) for each sample and cluster
df_avg_library_size <- colData(spe) %>%
    as.data.frame() %>%
    filter(PRECAST_cluster %in% wm_clusters) %>%
    group_by(sample_id, PRECAST_cluster) %>%
    summarise(avg_library_size = mean(sizeFactor))

# View the resulting dataframe
head(df_avg_library_size)

png(here("plots", "20_WM_comparisons", "sizeFactor_comparison.png"), unit="in", res=300, width = 7, height = 5)
ggplot(df_avg_library_size, aes(x = sample_id, y = avg_library_size, color = PRECAST_cluster)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("3"="blue","7"="brown","8"="grey"), name="PRECAST Cluster") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "Average Library Size by Sample and Cluster (WM regions)",
         x = "Sample ID",
         y = "Average Library Size")
dev.off()
