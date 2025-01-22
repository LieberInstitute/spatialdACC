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

dat_DLPFC <- as.data.frame(colData(sce_DLPFC))
dat_DLPFC <- dat_DLPFC[which(sce_DLPFC$layer_annotation %in% c("L5","L5/6")),]
summary_DLPFC <- dat_DLPFC %>%
    group_by(Sample) %>%
    summarize(total_DLPFC = n(), count38_DLPFC = sum(nmf38 > 0), count61_DLPFC = sum(nmf61 > 0))

summary_DLPFC$frac38_DLPFC <- summary_DLPFC$count38_DLPFC / summary_DLPFC$total_DLPFC
summary_DLPFC$frac61_DLPFC <- summary_DLPFC$count61_DLPFC / summary_DLPFC$total_DLPFC

dat_dACC <- as.data.frame(colData(sce))
dat_dACC <- dat_dACC[which(dat_dACC$cellType_azimuth %in% c("L5_IT")),]
summary_dACC <- dat_dACC %>%
    group_by(Sample) %>%
    summarize(total_dACC = n(), count38_dACC = sum(NMF_38 > 0), count61_dACC = sum(NMF_61 > 0))

summary_dACC$frac38_dACC <- summary_dACC$count38_dACC / summary_dACC$total_dACC
summary_dACC$frac61_dACC <- summary_dACC$count61_dACC / summary_dACC$total_dACC

dat_dACC <- as.data.frame(colData(sce))
dat_dACC <- dat_dACC[which(dat_dACC$cellType_azimuth %in% c("L5_ET")),]
summary_dACC <- dat_dACC %>%
    group_by(Sample) %>%
    summarize(total_dACC = n(), count38_dACC = sum(NMF_38 > 0), count61_dACC = sum(NMF_61 > 0))

summary_dACC$frac38_dACC <- summary_dACC$count38_dACC / summary_dACC$total_dACC
summary_dACC$frac61_dACC <- summary_dACC$count61_dACC / summary_dACC$total_dACC

plot_list <- list()

for (i in c(38,61)){
    print(paste0("i=", i))

    p <- plotColData(sce_DLPFC, x = "layer_annotation", y = paste0("nmf", i)) +
        ggtitle(paste0("NMF ", i, " DLPFC snRNA-seq Layer Boxplots")) +
        facet_wrap(~ sce_DLPFC$layer_annotation, scales = "free_x", nrow = 1) +
        ylim(c(0,0.0007))

    plot_list[[i]] <- p

}

pdf(file = here::here("plots", "13_NMF", "NMF_boxplots_DLPFC_single_nucleus_38_61.pdf"),
    width = 10, height = 10)

grid.arrange(
    grobs = plot_list[c(38,61)],
    ncol = 1,
    nrow = 2
    )

dev.off()
