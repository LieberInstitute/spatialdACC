setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(pheatmap)
library(scater)
library(scran)

# create violin plots to show the distribution of VAT1L expression in each cluster
load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))

sce <- computeSumFactors(sce)
sce <- logNormCounts(sce)

sce$counts_VAT1L <- counts(sce)[which(rowData(sce)$gene_name=="VAT1L"),]
sce$logcounts_VAT1L <- logcounts(sce)[which(rowData(sce)$gene_name=="VAT1L"),]

# create violin plots to show the distribution of VAT1L expression in each spatial domain
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata"))
spe$counts_VAT1L <- counts(spe)[which(rowData(spe)$gene_name=="VAT1L"),]
spe$logcounts_VAT1L <- logcounts(spe)[which(rowData(spe)$gene_name=="VAT1L"),]

pdf(file = here::here("plots", "16_VENs_analysis", "VAT1L_violin_plots.pdf"),
    width = 21, height = 20)

plotColData(sce, x = "cellType_azimuth", y = "counts_VAT1L") +
    ggtitle("VAT1L Counts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(spe, x = "layer", y = "counts_VAT1L") +
    ggtitle("VAT1L Counts by PRECAST Cluster") +
    facet_wrap(~ spe$layer, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "logcounts_VAT1L") +
    ggtitle("VAT1L Logcounts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(spe, x = "layer", y = "logcounts_VAT1L") +
    ggtitle("VAT1L Logcounts by PRECAST Cluster") +
    facet_wrap(~ spe$layer, scales = "free_x", nrow = 1)

dev.off()

# create violin plots to show the distribution of VAT1L expression in each cluster
sce_path_zip <- fetch_data("spatialDLPFC_snRNAseq")
sce_path <- unzip(sce_path_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
    file.path(tempdir(), "sce_DLPFC_annotated")
)

assay(sce, "logcounts") <- as(assay(sce, "logcounts"), "dgCMatrix")

sce$counts_VAT1L <- counts(sce)[which(rowData(sce)$gene_name=="VAT1L"),]
sce$logcounts_VAT1L <- logcounts(sce)[which(rowData(sce)$gene_name=="VAT1L"),]

# create violin plots to show the distribution of VAT1L expression in each spatial domain
spe_DLPFC_30 <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium")

spe_DLPFC_30$counts_VAT1L <- counts(spe_DLPFC_30)[which(rowData(spe_DLPFC_30)$gene_name=="VAT1L"),]
spe_DLPFC_30$logcounts_VAT1L <- logcounts(spe_DLPFC_30)[which(rowData(spe_DLPFC_30)$gene_name=="VAT1L"),]

# create spatial labels for DLPFC_30
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 3] <- "L2"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 8] <- "L4"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 7] <- "L6"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 5] <- "L3"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 6] <- "WM"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 4] <- "L5"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 2] <- "L1"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 1] <- "meninges"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 9] <- "WM"

spe_DLPFC_30$cluster <- spe_DLPFC_30$BayesSpace_harmony_09

sce_noNA <- sce[,!is.na(sce$layer_annotation)]

pdf(file = here::here("plots", "16_VENs_analysis", "VAT1L_violin_plots_DLPFC_30.pdf"),
    width = 21, height = 20)

plotColData(sce_noNA, x = "layer_annotation", y = "counts_VAT1L") +
    ggtitle("Single Nucleus VAT1L Counts by Layer Annotation") +
    facet_wrap(~ sce_noNA$layer_annotation, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_broad_hc", y = "counts_VAT1L") +
    ggtitle("Single Nucleus VAT1L Counts by Broad HC") +
    facet_wrap(~ sce$cellType_broad_hc, scales = "free_x", nrow = 1)

plotColData(spe_DLPFC_30, x = "cluster", y = "counts_VAT1L") +
    ggtitle("Visium VAT1L Counts by Cluster") +
    facet_wrap(~ spe_DLPFC_30$cluster, scales = "free_x", nrow = 1)

plotColData(sce_noNA, x = "layer_annotation", y = "logcounts_VAT1L") +
    ggtitle("Single Nucleus VAT1L Logcounts by Layer Annotation") +
    facet_wrap(~ sce_noNA$layer_annotation, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_broad_hc", y = "logcounts_VAT1L") +
    ggtitle("Single Nucleus VAT1L Logcounts by Broad HC") +
    facet_wrap(~ sce$cellType_broad_hc, scales = "free_x", nrow = 1)

plotColData(spe_DLPFC_30, x = "cluster", y = "logcounts_VAT1L") +
    ggtitle("Visium VAT1L Logcounts by Cluster") +
    facet_wrap(~ spe_DLPFC_30$cluster, scales = "free_x", nrow = 1)

dev.off()
