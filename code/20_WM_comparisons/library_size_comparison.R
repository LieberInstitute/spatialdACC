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
library(spatialLIBD)
library(patchwork)

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


p1 <- ggplot(df_avg_library_size, aes(x = sample_id, y = avg_library_size, color = PRECAST_cluster)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("3"="#377eb8","7"="#a65628","8"="#999999"), name="PRECAST Cluster") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "Average Library Size by Sample and Cluster (WM regions)",
         x = "Sample ID",
         y = "Average Library Size") +
    theme(legend.position="none")

load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))
colnames(spe) <- spe$key

precast_name <- paste0("nnSVG_PRECAST_captureArea_", 9)
spe <- cluster_import(
    spe,
    cluster_dir = here::here("processed-data", "08_clustering", "PRECAST", precast_name),
    prefix = paste0("clust", 9, "_")
    )

p2 <- vis_clus(spe = spe, sampleid = "V12N28-334_C1",
              clustervar = paste0("clust",9,"_PRECAST_cluster"), point_size = 0.9,
              spatial=F) +
    guides(fill = guide_legend(override.aes = list(size = 5))) +
    labs(fill="Cluster") +
    theme(legend.text = element_text(size = 10),legend.title = element_text(size=10)) +
    ggtitle("")

layout <- "
AAAB
AAAB
AAAB
AAAB
CCCB
"

png(here("plots", "20_WM_comparisons", "sizeFactor_comparison.png"), unit="in", res=300, width = 11, height = 8)
(p1 + p2) +
    plot_layout(design = layout)
dev.off()
