setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(pheatmap)
library(scater)
library(scran)
library(patchwork)

# create violin plots to show the distribution of RORB expression in each cluster
load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))

sce <- computeSumFactors(sce)
sce <- logNormCounts(sce)

sce$counts_RORB <- counts(sce)[which(rowData(sce)$gene_name=="RORB"),]
sce$logcounts_RORB <- logcounts(sce)[which(rowData(sce)$gene_name=="RORB"),]

sce_dACC <- sce

# create violin plots to show the distribution of RORB expression in each cluster
sce_path_zip <- fetch_data("spatialDLPFC_snRNAseq")
sce_path <- unzip(sce_path_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
    file.path(tempdir(), "sce_DLPFC_annotated")
)

assay(sce, "logcounts") <- as(assay(sce, "logcounts"), "dgCMatrix")

sce$counts_RORB <- counts(sce)[which(rowData(sce)$gene_name=="RORB"),]
sce$logcounts_RORB <- logcounts(sce)[which(rowData(sce)$gene_name=="RORB"),]

sce_dlPFC <- sce

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "celltype_colors.Rdata"))

p1 <- plotColData(sce_dACC, x = "cellType_azimuth", y = "logcounts_RORB", colour_by = "cellType_azimuth") +
    ggtitle("dACC RORB Logcounts by Azimuth Cell Type") +
    scale_colour_manual(values = celltype_colors) +
    facet_wrap(~ sce_dACC$cellType_azimuth, scales = "free_x", nrow = 1) +
    xlab("dACC Cell Type") +
    ylab("logcounts RORB") +
    theme(legend.position = "none")

p2 <- plotColData(sce_dlPFC, x = "cellType_layer", y = "logcounts_RORB") +
    ggtitle("dlPFC RORB Counts by Cell Type Layer") +
    facet_wrap(~ sce_dlPFC$cellType_layer, scales = "free_x", nrow = 1) +
    xlab("dlPFC Cell Type") +
    ylab("logcounts RORB") +
    theme(legend.position = "none")

png(file = here::here("plots", "11_differential_expression", "pseudobulk",
                      "boxplots_annotations", "RORB_boxplots_snRNA-seq.png"),
    width = 12, height = 12, unit="in", res = 300)

wrap_plots(p1,p2,ncol=1)

dev.off()
