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
library("svglite")
library("dplyr")
library("patchwork")
library("circlize")
library("Matrix")
library("dplyr")
library("tibble")

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))

# visualize heatmap of 2 genes for each Azimuth cell type
genes <- c("CUX2", "ADCYAP1", "VAT1L", "GABRQ", "RORB", "NTRK2", "PCP4","HTR2C",
           "FOXP2", "SEMA5A", "THEMIS", "CCN4", "CUX1", "SEMA6D", "ZFHX3", "NXPH4")

cell_type_per_gene <- c("L2_3_IT", "L2_3_IT", "L5_ET", "L5_ET", "L5_IT", "L5_IT",
                        "L5_6_NP", "L5_6_NP", "L6_CT", "L6_CT", "L6_IT", "L6_IT",
                        "L6_IT_Car3", "L6_IT_Car3", "L6b", "L6b")

# only keep cell types that are in the genes list
sce <- sce[, sce$cellType_azimuth %in% cell_type_per_gene]
sce$cellType_azimuth <- as.factor(sce$cellType_azimuth)

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "celltype_colors.Rdata"))

sce <- logNormCounts(sce)

log_mat <- logcounts(sce)[genes, , drop = FALSE]
cell_types <- colData(sce)$cellType_azimuth

long_expr <- as.data.frame(as.matrix(log_mat)) %>%
    rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "cell", values_to = "logcounts") %>%
    mutate(cell_type = cell_types[match(cell, colnames(sce))])

mean_expr <- long_expr %>%
    group_by(gene, cell_type) %>%
    summarise(mean_logcounts = mean(logcounts), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = cell_type, values_from = mean_logcounts)

heatmap_mat <- as.matrix(column_to_rownames(mean_expr, var = "gene"))

# reorder based on genes variable
heatmap_mat <- heatmap_mat[genes, ]
# make second column fourth column
heatmap_mat <- heatmap_mat[, c(1, 3, 4, 2, 5:ncol(heatmap_mat))]

heatmap_mat <- t(scale(t(heatmap_mat)))

pdf(here("plots", "snRNA-seq", "05_azimuth", "azimuth_top_genes_heatmap.pdf"), height = 3, width = 5)

Heatmap(
    t(heatmap_mat),
    cluster_rows = F, cluster_columns = F, row_names_side = "left",
    show_row_names = TRUE, show_column_names = TRUE,
    heatmap_legend_param = list(
        title = "scaled avg.\nlogcount\nexpression", at = c(-4, 0, 4),
        labels = c("-4", "0", "4")
    ),
    top_annotation = HeatmapAnnotation(
        " " = cell_type_per_gene,
        col = list(" " = celltype_colors),
        show_legend = F
    ),
)

dev.off()
