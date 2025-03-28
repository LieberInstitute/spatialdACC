setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library(here)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(gridExtra)
library(grid)

load(file=here("processed-data","19_bulk_deconvolution","combined_pvalues_ordered_spatial_pain.Rdata"))
load(file=here("processed-data","19_bulk_deconvolution","combined_pvalues_ordered_sn_pain.Rdata"))
load(file=here("processed-data","19_bulk_deconvolution","combined_pvalues_ordered_spatial_bulk.Rdata"))
load(file=here("processed-data","19_bulk_deconvolution","combined_pvalues_ordered_sn_bulk.Rdata"))

row.names(combined_pvalues_ordered_spatial_bulk)[1] <- "           L1"
row.names(combined_pvalues_ordered_spatial_pain)[1] <- "          L6b"
colnames(combined_pvalues_ordered_spatial_pain)[1] <- "Pain Down"
colnames(combined_pvalues_ordered_sn_pain)[1] <- "Pain Down"
colnames(combined_pvalues_ordered_spatial_pain)[2] <- "Pain Up"
colnames(combined_pvalues_ordered_sn_pain)[2] <- "Pain Up"

common_height <- unit(5, "mm")
common_width <- unit(5, "mm")

col_fun_1 <- colorRamp2(
    c(1.3, max(-log10(combined_pvalues_ordered_sn_bulk))),
    c("white", "blue") # White for -log10(p) >= 1.3 (p >= 0.05), blue for more significant p-values
)

p1 <- Heatmap(
    -log10(combined_pvalues_ordered_sn_bulk),
    name = "-log(p)",
    col = col_fun_1,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_side = "left",
    column_title = "",
    heatmap_legend_param = list(
        title = "-log(p)",
        title_position = "topcenter",
        title_gp = gpar(fontsize = 10),
        labels_gp = gpar(fontsize = 8)
    ),
    column_names_side = "bottom",
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 2),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    height = common_height * nrow(combined_pvalues_ordered_sn_bulk),
    width = common_width * ncol(combined_pvalues_ordered_sn_bulk)
)


col_fun_2 <- colorRamp2(
    c(1.3, max(-log10(combined_pvalues_ordered_spatial_bulk))),
    c("white", "blue") # White for -log10(p) >= 1.3 (p >= 0.05), blue for more significant p-values
)


p2 <- Heatmap(
    -log10(combined_pvalues_ordered_spatial_bulk),
    name = "-log(p)",
    col = col_fun_2,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_side = "left",
    column_title = "",
    heatmap_legend_param = list(
        title = "-log(p)",
        title_position = "topcenter",
        title_gp = gpar(fontsize = 10),
        labels_gp = gpar(fontsize = 8)
    ),
    column_names_side = "bottom",
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 2),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    height = common_height * nrow(combined_pvalues_ordered_spatial_bulk),
    width = common_width * ncol(combined_pvalues_ordered_spatial_bulk)
)

col_fun_3 <- colorRamp2(
    c(1.3, max(-log10(combined_pvalues_ordered_sn_pain))),
    c("white", "blue") # White for -log10(p) >= 1.3 (p >= 0.05), blue for more significant p-values
)

p3 <- Heatmap(
    -log10(combined_pvalues_ordered_sn_pain),
    name = "-log(p)",
    col = col_fun_3,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_side = "left",
    column_title = "",
    heatmap_legend_param = list(
        title = "-log(p)",
        title_position = "topcenter",
        title_gp = gpar(fontsize = 10),
        labels_gp = gpar(fontsize = 8)
    ),
    column_names_side = "bottom",
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 2),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    height = common_height * nrow(combined_pvalues_ordered_sn_pain),
    width = common_width * ncol(combined_pvalues_ordered_sn_pain)
)


col_fun_4 <- colorRamp2(
    c(1.3, max(-log10(combined_pvalues_ordered_spatial_pain))),
    c("white", "blue") # White for -log10(p) >= 1.3 (p >= 0.05), blue for more significant p-values
)


p4 <- Heatmap(
    -log10(combined_pvalues_ordered_spatial_pain),
    name = "-log(p)",
    col = col_fun_4,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_side = "left",
    column_title = "",
    heatmap_legend_param = list(
        title = "-log(p)",
        title_position = "topcenter",
        title_gp = gpar(fontsize = 10),
        labels_gp = gpar(fontsize = 8)
    ),
    column_names_side = "bottom",
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 2),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    height = common_height * nrow(combined_pvalues_ordered_spatial_pain),
    width = common_width * ncol(combined_pvalues_ordered_spatial_pain)
)

pdf(here("plots", "21_pain_enrichment", "enrichment_heatmaps.pdf"), height = 8.2, width = 4.8)

grid.newpage()
pushViewport(viewport(layout = grid.layout(
    3, 2,
    heights = unit(c(1, 7, 19), "null"),
    widths = unit(c(4, 2), "null")
)))

define_position <- function(plot, row, col) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    draw(plot, newpage = FALSE)
    popViewport()
}

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
grid.text("Clinical and Pain Enrichment", gp = gpar(fontsize = 16, fontface = "bold"))
popViewport()

define_position(p2, 2, 1)
define_position(p4, 2, 2)
define_position(p1, 3, 1)
define_position(p3, 3, 2)

dev.off()
