setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("SingleCellExperiment")
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

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))

#plot harmony UMAP colored by azimuth labels
pdf(file = here::here("plots", "snRNA-seq", "05_azimuth", "HARMONY_azimuth_UMAP.pdf"))

sce$cellType_azimuth <- as.factor(sce$cellType_azimuth)
plotReducedDim(sce, dimred="UMAP-HARMONY", colour_by="cellType_azimuth", point_size = 0.5)

dev.off()


