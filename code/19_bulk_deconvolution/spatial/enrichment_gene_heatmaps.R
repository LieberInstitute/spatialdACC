setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library(here)
library(ComplexHeatmap)
library(circlize)
library(EnhancedVolcano)
library(SpatialExperiment)
library(RColorBrewer)

# Load DEG data
load(
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_DE", "nnSVG_PRECAST_captureArea_9.Rdata")
)

nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", 9)
load(
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_pseudobulk", paste0(nnSVG_precast_name,".Rdata"))
)

# Subset modeling results
enrichment_results <- modeling_results[["enrichment"]]

# Load the bulk data
file = here("processed-data", "PTSD_bulk", "appi.ajp.21020162.ds003.csv")
bulk <- read.csv(file, row.names = 1)

# Step 1: Extract base gene IDs (remove version numbers)
gene_ids <- rownames(bulk)
base_gene_ids <- sub("\\..*", "", gene_ids)

# Step 2: Check for duplicates
duplicated_genes <- base_gene_ids[duplicated(base_gene_ids)]

# Print duplicated genes (if any)
if (length(duplicated_genes) > 0) {
    print("Duplicated genes found (without version numbers):")
    print(duplicated_genes)
} else {
    print("No duplicated genes found (without version numbers).")
}

# Remove .xx from gene names
rownames(bulk) <- base_gene_ids

# List columns of interest from bulk data
vars <- c("dACC_logFC_PTSD", "dACC_t_PTSD", "dACC_adjPVal_PTSD",
          "dACC_logFC_MDD", "dACC_t_MDD", "dACC_adjPVal_MDD")

# Create smaller df with these vars
bulk <- bulk[, vars]

# Overlap to get genes in both datasets
overlap <- intersect(rownames(bulk), rownames(enrichment_results))
length(overlap) # [1] 12203

# Subset the bulk object
bulk <- bulk[overlap, ]

# Subset the modeling results
enrichment_results <- enrichment_results[overlap, ]

# add the gene names to the bulk dataset because they are already in the same order
bulk$gene_name <- enrichment_results$gene

top_n <- 500
plot_list_500 <- list()

for (i in c("L1","L2","L3","L5","L6a","L6b","WM")) {
    # choose the top n by t statistics
    t_stat_threshold <- sort(enrichment_results[[paste0("t_stat_", i)]], decreasing = T)[top_n]
    DE_clust_genes_up <- rownames(enrichment_results[enrichment_results[[paste0("t_stat_", i)]] >= t_stat_threshold, ])

    # make heatmap for this layer's markers only
    spe_sub <- spe_pseudo[which(rownames(spe_pseudo) %in% DE_clust_genes_up),]
    expr_mat <- logcounts(spe_sub)
    expr_mat <- t(scale(t(expr_mat)))
    qs <- quantile(expr_mat, c(0.1, 0.95))
    col_fun <- circlize::colorRamp2(c(qs[1], 0, qs[2]), c("#FF00FF", "black", "#FFFF00"))
    cluster_anno <- spe_sub$layer

    p1 <- Heatmap(expr_mat, name = paste0("Top 500 Markers ",i),
                  column_split = factor(cluster_anno),
                  cluster_columns = TRUE,
                  show_column_dend = FALSE,
                  cluster_column_slices = TRUE,
                  column_title_gp = gpar(fontsize = 8),
                  column_gap = unit(0.5, "mm"),
                  cluster_rows = TRUE,
                  show_row_dend = FALSE,
                  col = col_fun,
                  row_names_gp = gpar(fontsize = 4),
                  column_title_rot = 90,
                  top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
                  show_column_names = FALSE,
                  use_raster = TRUE,
                  raster_quality = 4)


    plot_list_500[[i]] <- p1
}

pdf(file=here("plots","19_bulk_deconvolution","spatial_genes_top_500.pdf"))
print(plot_list_500[["L1"]])
print(plot_list_500[["L2"]])
print(plot_list_500[["L3"]])
print(plot_list_500[["L5"]])
print(plot_list_500[["L6a"]])
print(plot_list_500[["L6b"]])
print(plot_list_500[["WM"]])
dev.off()

plot_list_adjp <- list()

for (i in c("L1","L2","L3","L5","L6a","L6b","WM")) {

    DE_clust_genes_up <- rownames(enrichment_results[
        enrichment_results[[paste0("fdr_", i)]] < 0.05 & enrichment_results[[paste0("logFC_", i)]] > 1.5, ])

    # make heatmap for this layer's markers only
    spe_sub <- spe_pseudo[which(rownames(spe_pseudo) %in% DE_clust_genes_up),]
    expr_mat <- logcounts(spe_sub)
    expr_mat <- t(scale(t(expr_mat)))
    qs <- quantile(expr_mat, c(0.1, 0.95))
    col_fun <- circlize::colorRamp2(c(qs[1], 0, qs[2]), c("#FF00FF", "black", "#FFFF00"))
    cluster_anno <- spe_sub$layer

    p1 <- Heatmap(expr_mat, name = paste0("adj p < 0.05 & LFC > 1.5 ",i),
                  column_split = factor(cluster_anno),
                  cluster_columns = TRUE,
                  show_column_dend = FALSE,
                  cluster_column_slices = TRUE,
                  column_title_gp = gpar(fontsize = 8),
                  column_gap = unit(0.5, "mm"),
                  cluster_rows = TRUE,
                  show_row_dend = FALSE,
                  col = col_fun,
                  row_names_gp = gpar(fontsize = 4),
                  column_title_rot = 90,
                  top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
                  show_column_names = FALSE,
                  use_raster = TRUE,
                  raster_quality = 4)

    plot_list_adjp[[i]] <- p1

}

pdf(file=here("plots","19_bulk_deconvolution","spatial_genes_0.05_LFC_1.5.pdf"))
print(plot_list_adjp[["L1"]])
print(plot_list_adjp[["L2"]])
print(plot_list_adjp[["L3"]])
print(plot_list_adjp[["L5"]])
print(plot_list_adjp[["L6a"]])
print(plot_list_adjp[["L6b"]])
print(plot_list_adjp[["WM"]])
dev.off()
