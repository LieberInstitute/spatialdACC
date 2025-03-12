setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library(here)
library(ComplexHeatmap)
library(circlize)
library(EnhancedVolcano)
library(SpatialExperiment)
library(RColorBrewer)

# load DEGs
load(
    file = here("processed-data", "snRNA-seq", "05_azimuth", "azimuth_DE_results.Rdata")
)

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_pseudo_azimuth_pseudo.Rdata"))

# replace spaces with underscores
colData(sce_pseudo)$cellType_azimuth <- replace(colData(sce_pseudo)$cellType_azimuth, colData(sce_pseudo)$cellType_azimuth == "L2/3 IT", "L2_3_IT")
colData(sce_pseudo)$cellType_azimuth <- replace(colData(sce_pseudo)$cellType_azimuth, colData(sce_pseudo)$cellType_azimuth == "L5 ET", "L5_ET")
colData(sce_pseudo)$cellType_azimuth <- replace(colData(sce_pseudo)$cellType_azimuth, colData(sce_pseudo)$cellType_azimuth == "L5 IT", "L5_IT")
colData(sce_pseudo)$cellType_azimuth <- replace(colData(sce_pseudo)$cellType_azimuth, colData(sce_pseudo)$cellType_azimuth == "L5/6 NP", "L5_6_NP")
colData(sce_pseudo)$cellType_azimuth <- replace(colData(sce_pseudo)$cellType_azimuth, colData(sce_pseudo)$cellType_azimuth == "L6 CT", "L6_CT")
colData(sce_pseudo)$cellType_azimuth <- replace(colData(sce_pseudo)$cellType_azimuth, colData(sce_pseudo)$cellType_azimuth == "L6 IT", "L6_IT")
colData(sce_pseudo)$cellType_azimuth <- replace(colData(sce_pseudo)$cellType_azimuth, colData(sce_pseudo)$cellType_azimuth == "L6 IT Car3", "L6_IT_Car3")
colData(sce_pseudo)$cellType_azimuth <- replace(colData(sce_pseudo)$cellType_azimuth, colData(sce_pseudo)$cellType_azimuth == "Sst Chodl", "Sst_Chodl")
colData(sce_pseudo)$cellType_azimuth <- replace(colData(sce_pseudo)$cellType_azimuth, colData(sce_pseudo)$cellType_azimuth == "MicroPVM", "Micro_PVM")

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

k <- unique(sce_pseudo$cellType_azimuth)

top_n <- 500
plot_list_500 <- list()

for (i in k) {
    # choose the top n by t statistics
    t_stat_threshold <- sort(enrichment_results[[paste0("t_stat_", i)]], decreasing = T)[top_n]
    DE_clust_genes_up <- rownames(enrichment_results[enrichment_results[[paste0("t_stat_", i)]] >= t_stat_threshold, ])

    # make heatmap for this layer's markers only
    sce_sub <- sce_pseudo[which(rowData(sce_pseudo)$gene_id %in% DE_clust_genes_up),]
    expr_mat <- logcounts(sce_sub)
    expr_mat <- t(scale(t(expr_mat)))
    qs <- quantile(expr_mat, c(0.1, 0.95))
    col_fun <- circlize::colorRamp2(c(qs[1], 0, qs[2]), c("#FF00FF", "black", "#FFFF00"))
    cluster_anno <- sce_sub$cellType_azimuth

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

pdf(file=here("plots","19_bulk_deconvolution","single_nucleus_genes_top_500.pdf"))
print(plot_list_500[["Astro"]])
print(plot_list_500[["Endo"]])
print(plot_list_500[["L2_3_IT"]])
print(plot_list_500[["L5_6_NP"]])
print(plot_list_500[["L5_ET"]])
print(plot_list_500[["L5_IT"]])
print(plot_list_500[["L6_CT"]])
print(plot_list_500[["L6_IT"]])
print(plot_list_500[["L6_IT_Car3"]])
print(plot_list_500[["L6b"]])
print(plot_list_500[["Lamp5"]])
print(plot_list_500[["Micro_PVM"]])
print(plot_list_500[["Oligo"]])
print(plot_list_500[["OPC"]])
print(plot_list_500[["Pvalb"]])
print(plot_list_500[["Sncg"]])
print(plot_list_500[["Sst"]])
print(plot_list_500[["Vip"]])
print(plot_list_500[["VLMC"]])

dev.off()

plot_list_adjp <- list()

for (i in k) {

    DE_clust_genes_up <- rownames(enrichment_results[
        enrichment_results[[paste0("fdr_", i)]] < 0.05 & enrichment_results[[paste0("logFC_", i)]] > 2, ])

    # make heatmap for this layer's markers only
    sce_sub <- sce_pseudo[which(rowData(sce_pseudo)$gene_id %in% DE_clust_genes_up),]
    expr_mat <- logcounts(sce_sub)
    expr_mat <- t(scale(t(expr_mat)))
    qs <- quantile(expr_mat, c(0.1, 0.95))
    col_fun <- circlize::colorRamp2(c(qs[1], 0, qs[2]), c("#FF00FF", "black", "#FFFF00"))
    cluster_anno <- sce_sub$cellType_azimuth

    p1 <- Heatmap(expr_mat, name = paste0("adj p < 0.05 & LFC > 2 ",i),
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

pdf(file=here("plots","19_bulk_deconvolution","single_nucleus_genes_0.05_LFC_2.pdf"))
print(plot_list_adjp[["Astro"]])
print(plot_list_adjp[["Endo"]])
print(plot_list_adjp[["L2_3_IT"]])
print(plot_list_adjp[["L5_6_NP"]])
print(plot_list_adjp[["L5_ET"]])
print(plot_list_adjp[["L5_IT"]])
print(plot_list_adjp[["L6_CT"]])
print(plot_list_adjp[["L6_IT"]])
print(plot_list_adjp[["L6_IT_Car3"]])
print(plot_list_adjp[["L6b"]])
print(plot_list_adjp[["Lamp5"]])
print(plot_list_adjp[["Micro_PVM"]])
print(plot_list_adjp[["Oligo"]])
print(plot_list_adjp[["OPC"]])
print(plot_list_adjp[["Pvalb"]])
print(plot_list_adjp[["Sncg"]])
print(plot_list_adjp[["Sst"]])
print(plot_list_adjp[["Vip"]])
print(plot_list_adjp[["VLMC"]])
dev.off()
