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

load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata"))

# compare library size of WM1, WM2, WM-CC
# overall and separately by sample
# make a dot plot of the library size of each sample

# create a dataframe by taking the average of colData(spe)$sizeFactor for each sample in colData(spe)$sample_id and cluster in colData(spe)$PRECAST_cluster

# Filter for relevant clusters
wm_clusters <- c("WM1", "WM2", "WM-CC")

# Create a dataframe with the average library size (sizeFactor) for each sample and cluster
df_avg_library_size <- colData(spe) %>%
    as.data.frame() %>%
    filter(PRECAST_cluster %in% wm_clusters) %>%
    group_by(sample_id, PRECAST_cluster) %>%
    summarise(avg_library_size = mean(sizeFactor))

# View the resulting dataframe
head(df_avg_library_size)

pdf(here("plots", "20_WM_comparisons", "sizeFactor_comparison.pdf"), width = 20, height = 5)
ggplot(df_avg_library_size, aes(x = sample_id, y = avg_library_size, color = PRECAST_cluster)) +
    geom_point(size = 3) +
    # make x axis labels vertical
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme_bw() +
    labs(title = "Average Library Size by Sample and Cluster (WM1, WM2, WM-CC)",
         x = "Sample ID",
         y = "Average Library Size")
dev.off()
