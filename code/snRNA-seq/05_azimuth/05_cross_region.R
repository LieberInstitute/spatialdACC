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

# use cellType_layer annotation
sum_by_sample <- setNames(aggregate(sum ~ Sample, colData(sce), sum), c("Sample", "sum_sample"))
detected_by_sample <- setNames(aggregate(detected ~ Sample, colData(sce), sum), c("Sample", "detected_sample"))

sce_pseudo <-
    registration_pseudobulk(sce,
                            var_registration = "cellType_layer",
                            var_sample_id = "Sample"
    )

pca <- prcomp(t(assays(sce_pseudo)$logcounts))
pca_pseudo <- pca$x[, seq_len(20)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(sce_pseudo) <- list(PCA = pca_pseudo)

## Plot PCs
col_data_df <- as.data.frame(colData(sce_pseudo))
col_data_df <- left_join(col_data_df, detected_by_sample, by = "Sample")
col_data_df <- left_join(col_data_df, sum_by_sample, by = "Sample")

colData(sce_pseudo) <- DataFrame(col_data_df)

# make supp figure
p1 <- plotPCA(
    sce_pseudo,
    colour_by = "cellType_layer",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(sce_pseudo)$PCA_var_explained
)

p2 <- plotPCA(
    sce_pseudo,
    colour_by = "Sample",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(sce_pseudo)$PCA_var_explained
)

p3 <- plotPCA(
    sce_pseudo,
    colour_by = "sum_sample",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(sce_pseudo)$PCA_var_explained
)

vars <- getVarianceExplained(sce_pseudo,
                             variables = c("cellType_layer","Sample", "sum_sample", "detected_sample")
)

p4 <- plotExplanatoryVariables(vars)

png(file = here("plots", "snRNA-seq",
                "05_azimuth",
                "pseudobulk_PC_DLPFC.png"),
    width = 10, height = 10, unit="in", res=300)

wrap_plots(p1,p2,p3,p4,nrow=2)
dev.off()

save(sce_pseudo, file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_dlPFC_pseudo.Rdata"))

modeling_results <- registration_wrapper(
    sce,
    var_registration = "cellType_layer",
    var_sample_id = "Sample",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)


# save modeling results list
save(
    modeling_results,
    file = here("processed-data", "snRNA-seq", "05_azimuth", "dlPFC_DE_results.Rdata")
)

# spatial registration with Azimuth cell types in dACC dataset
## extract t-statics and rename
registration_t_stats <- modeling_results$enrichment[, grep("^t_stat", colnames(modeling_results$enrichment))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))

dim(registration_t_stats)

## check out table
registration_t_stats[1:5, 1:5]

# load Azimuth modeling results list
load(
    file = here("processed-data", "snRNA-seq", "05_azimuth", "azimuth_DE_results.Rdata")
)

cor_layer <- layer_stat_cor(
    stats = registration_t_stats,
    modeling_results = modeling_results,
    model_type = "enrichment",
    top_n = 100
)

anno <- annotate_registered_clusters(
    cor_stats_layer = cor_layer,
    confidence_threshold = 0.25,
    cutoff_merge_ratio = 0.25
)

anno

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "celltype_colors.Rdata"))


pdf(file = here::here("plots", "12_spatial_registration", "azimuth",
                      paste0("azimuth_","dlPFC","_heatmap.pdf")), width = 6, height = 6)
layer_stat_cor_plot(cor_layer, color_max = max(cor_layer),
                    reference_colors = celltype_colors,
                    annotation = anno,
                    cluster_rows = FALSE,
                    cluster_columns = FALSE)
dev.off()

