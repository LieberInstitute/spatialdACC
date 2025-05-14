setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(scater)
library(RcppML)
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

# get NMF results from single nucleus data
x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))

sce_path_zip <- fetch_data("spatialDLPFC_snRNAseq")
sce_path <- unzip(sce_path_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
    file.path(tempdir(), "sce_DLPFC_annotated")
)
assay(sce, "logcounts") <- as(assay(sce, "logcounts"), "dgCMatrix")
sce_DLPFC <- sce

# load Single Nucleus object
load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))

# extract patterns
patterns <- t(x@h)
colnames(patterns) <- paste("NMF", 1:75, sep = "_")

loadings <- x@w
rownames(loadings) <- rownames(sce)

# ====== project loadings to spatial data =======
i <- intersect(rowData(sce_DLPFC)$gene_name,rownames(loadings))
loadings <- loadings[rownames(loadings) %in% i,]
sce_DLPFC <- sce_DLPFC[rowData(sce_DLPFC)$gene_name %in% i,]
loadings <- loadings[match(rowData(sce_DLPFC)$gene_name,rownames(loadings)),]

logcounts <- logcounts(sce_DLPFC)
#loadings <- as(loadings, "dgCMatrix")

proj <- project(w=loadings, data=logcounts)

# remove rowSums == 0
proj <- proj[rowSums(proj) != 0,]

proj <- apply(proj,1,function(x){x/sum(x)})

colData(sce_DLPFC) <- cbind(colData(sce_DLPFC),proj)
colData(sce) <- cbind(colData(sce),patterns)

# check if NMF38 and 61 were removed
# no

# save sce_DLPFC
save(sce_DLPFC, file = here("processed-data", "13_NMF", "sce_NMF_DLPFC_30_snRNA-seq.Rdata"))

sce_dACC <- sce

dat_DLPFC <- as.data.frame(colData(sce_DLPFC))
dat_dACC <- as.data.frame(colData(sce_dACC))

# save for future use
save(dat_dACC, dat_DLPFC, file = here("processed-data", "13_NMF", "DLPFC_dACC_celltype_NMF.Rdata"))


plot_list_DLPFC <- list()
sce_DLPFC <- sce_DLPFC[which(sce_DLPFC$cellType_azimuth %in% c("L5 IT", "L5 ET")),]


for (i in c(38,61)){
    print(paste0("i=", i))

    p <- plotColData(sce_DLPFC, x = "cellType_azimuth", y = paste0("nmf", i)) +
        ggtitle(paste0("NMF ", i, " DLPFC snRNA-seq Layer Boxplots")) +
        facet_wrap(~ sce_DLPFC$cellType_azimuth, scales = "free_x", nrow = 1) +
        ylim(c(0,0.004))

    plot_list_DLPFC[[i]] <- p

}

plot_list_dACC <- list()

for (i in c(38,61)){
    print(paste0("i=", i))

    p <- plotColData(sce_dACC, x = "cellType_azimuth", y = paste0("NMF_", i)) +
        ggtitle(paste0("NMF ", i, " DLPFC snRNA-seq Layer Boxplots")) +
        facet_wrap(~ sce_dACC$cellType_azimuth, scales = "free_x", nrow = 1) +
        ylim(c(0,0.004))

    plot_list_dACC[[i]] <- p

}

pdf(file = here::here("plots", "13_NMF", "NMF_boxplots_DLPFC_single_nucleus_38_61.pdf"),
    width = 10, height = 10)

grid.arrange(
    grobs = plot_list_DLPFC[c(38,61)],
    ncol = 1,
    nrow = 2
    )

grid.arrange(
    grobs = plot_list_dACC[c(38,61)],
    ncol = 1,
    nrow = 2
)

dev.off()
