setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(pheatmap)
library(scater)

# create violin plots to show the distribution of VAT1L expression in each cluster
load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))

sce <- computeSumFactors(sce)
sce <- logNormCounts(sce)

sce$counts_VAT1L <- counts(sce)[which(rowData(sce)$gene_name=="VAT1L"),]
sce$logcounts_VAT1L <- logcounts(sce)[which(rowData(sce)$gene_name=="VAT1L"),]


# create violin plots to show the distribution of VAT1L expression in each spatial domain
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9.Rdata"))
spe$counts_VAT1L <- counts(spe)[which(rowData(spe)$gene_name=="VAT1L"),]
spe$logcounts_VAT1L <- logcounts(spe)[which(rowData(spe)$gene_name=="VAT1L"),]

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

pdf(file = here::here("plots", "16_VENs_analysis", "VAT1L_violin_plots.pdf"),
    width = 21, height = 20)

plotColData(sce, x = "cellType_azimuth", y = "counts_VAT1L") +
    ggtitle("VAT1L Counts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(spe, x = "PRECAST_cluster", y = "counts_VAT1L") +
    ggtitle("VAT1L Counts by PRECAST Cluster") +
    facet_wrap(~ spe$PRECAST_cluster, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "logcounts_VAT1L") +
    ggtitle("VAT1L Logcounts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(spe, x = "PRECAST_cluster", y = "logcounts_VAT1L") +
    ggtitle("VAT1L Logcounts by PRECAST Cluster") +
    facet_wrap(~ spe$PRECAST_cluster, scales = "free_x", nrow = 1)

dev.off()
