setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(pheatmap)
library(scater)

# create violin plots to show the distribution of DRD5 expression in each cluster
load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))
sce$counts_DRD5 <- counts(sce)[which(rowData(sce)$gene_name=="DRD5"),]


# create violin plots to show the distribution of DRD5 expression in each spatial domain
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9.Rdata"))
spe$counts_DRD5 <- counts(spe)[which(rowData(spe)$gene_name=="DRD5"),]

spe$PRECAST_cluster <- unfactor(spe$PRECAST_cluster)
spe$PRECAST_cluster[spe$PRECAST_cluster == 3] <- "WM1"
spe$PRECAST_cluster[spe$PRECAST_cluster == 8] <- "WM2"
spe$PRECAST_cluster[spe$PRECAST_cluster == 7] <- "WM-CC"
spe$PRECAST_cluster[spe$PRECAST_cluster == 5] <- "L6b"
spe$PRECAST_cluster[spe$PRECAST_cluster == 6] <- "L6a"
spe$PRECAST_cluster[spe$PRECAST_cluster == 4] <- "L5"
spe$PRECAST_cluster[spe$PRECAST_cluster == 2] <- "L3"
spe$PRECAST_cluster[spe$PRECAST_cluster == 1] <- "L2"
spe$PRECAST_cluster[spe$PRECAST_cluster == 9] <- "L1"

pdf(file = here::here("plots", "16_VENs_analysis", "DRD5_violin_plots.pdf"),
    width = 21, height = 20)

plotColData(sce, x = "cellType_azimuth", y = "counts_DRD5") +
    ggtitle("DRD5 Counts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(spe, x = "PRECAST_cluster", y = "counts_DRD5") +
    ggtitle("DRD5 Counts by PRECAST Cluster") +
    facet_wrap(~ spe$PRECAST_cluster, scales = "free_x", nrow = 1)

dev.off()
