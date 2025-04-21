setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(scater)
library(ggspavis)
library(here)
library(scRNAseq)
library(Matrix)
library(scran)
library(scuttle)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(igraph)
library(bluster)
library(patchwork)
library(cowplot)
library(projectR)
library(spatialLIBD)
library(gridExtra)
library(dplyr)

sce_path_zip <- fetch_data("spatialDLPFC_snRNAseq")
sce_path <- unzip(sce_path_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
    file.path(tempdir(), "sce_DLPFC_annotated")
)
assay(sce, "logcounts") <- as(assay(sce, "logcounts"), "dgCMatrix")

# "Ambig" out of sce$cellType_broad_hc
idx <- which(sce$cellType_broad_hc == "Ambiguous")
sce <- sce[,-idx]
sce_dlPFC <- sce

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))
sce_dACC <- sce

GABRQ_dACC_SULF2_dACC_list <- list()
GABRQ_dACC_POU3F1_dACC_list <- list()
SULF2_dACC_POU3F1_dACC_list <- list()

GABRQ_dlPFC_SULF2_dlPFC_list <- list()
GABRQ_dlPFC_POU3F1_dlPFC_list <- list()
SULF2_dlPFC_POU3F1_dlPFC_list <- list()

for (br in unique(colData(sce_dACC)$brain)) {
    print(br)
    sce_dACC_subset <- sce_dACC[,which(colData(sce_dACC)$brain == br)]
    GABRQ_dACC <- counts(sce_dACC_subset)[which(rownames(counts(sce_dACC_subset))=="GABRQ"),]
    SULF2_dACC <- counts(sce_dACC_subset)[which(rownames(counts(sce_dACC_subset))=="SULF2"),]
    POU3F1_dACC <- counts(sce_dACC_subset)[which(rownames(counts(sce_dACC_subset))=="POU3F1"),]

    GABRQ_dACC_SULF2_dACC_list <- append(GABRQ_dACC_SULF2_dACC_list, cor(GABRQ_dACC, SULF2_dACC))
    GABRQ_dACC_POU3F1_dACC_list <- append(GABRQ_dACC_POU3F1_dACC_list, cor(GABRQ_dACC, POU3F1_dACC))
    SULF2_dACC_POU3F1_dACC_list <- append(SULF2_dACC_POU3F1_dACC_list, cor(SULF2_dACC, POU3F1_dACC))

    sce_dlPFC_subset <- sce_dlPFC[,which(colData(sce_dlPFC)$BrNum == br)]
    GABRQ_dlPFC <- counts(sce_dlPFC_subset)[which(rownames(counts(sce_dlPFC_subset))=="GABRQ"),]
    SULF2_dlPFC <- counts(sce_dlPFC_subset)[which(rownames(counts(sce_dlPFC_subset))=="SULF2"),]
    POU3F1_dlPFC <- counts(sce_dlPFC_subset)[which(rownames(counts(sce_dlPFC_subset))=="POU3F1"),]

    GABRQ_dlPFC_SULF2_dlPFC_list <- append(GABRQ_dlPFC_SULF2_dlPFC_list, cor(GABRQ_dlPFC, SULF2_dlPFC))
    GABRQ_dlPFC_POU3F1_dlPFC_list <- append(GABRQ_dlPFC_POU3F1_dlPFC_list, cor(GABRQ_dlPFC, POU3F1_dlPFC))
    SULF2_dlPFC_POU3F1_dlPFC_list <- append(SULF2_dlPFC_POU3F1_dlPFC_list, cor(SULF2_dlPFC, POU3F1_dlPFC))

}

corr_list <- c(unlist(GABRQ_dACC_SULF2_dACC_list),
               unlist(GABRQ_dACC_POU3F1_dACC_list),
               unlist(SULF2_dACC_POU3F1_dACC_list),
               unlist(GABRQ_dlPFC_SULF2_dlPFC_list),
               unlist(GABRQ_dlPFC_POU3F1_dlPFC_list),
               unlist(SULF2_dlPFC_POU3F1_dlPFC_list))

region_list <- c(rep("dACC",30), rep("dlPFC",30))

genes_list <- c(rep(c("GABRQ & SULF2"),10),
                rep(c("POU3F1 & GABRQ"),10),
                rep(c("POU3F1 & SULF2"),10),
                rep(c("GABRQ & SULF2"),10),
                rep(c("POU3F1 & GABRQ"),10),
                rep(c("POU3F1 & SULF2"),10))

df <- data.frame(
    Correlation = corr_list,
    Region = region_list,
    Genes = genes_list
)

p <- ggplot(df, aes(x=Genes, y=Correlation, color=Region)) +
    geom_point() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1))

png(here("plots", "snRNA-seq", "05_azimuth", "VENs_correlation_dotplot.png"), height=5, width=5, unit="in",res=300)
print(p)
dev.off()
