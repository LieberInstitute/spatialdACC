setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library(here)
library(ComplexHeatmap)
library(circlize)
library(SingleCellExperiment)
library(grid)

# load DEGs
load(
    file = here("processed-data", "snRNA-seq", "05_azimuth", "azimuth_DE_results.Rdata")
)

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))

# replace spaces with underscores
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L2/3 IT", "L2_3_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5 ET", "L5_ET")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5 IT", "L5_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5/6 NP", "L5_6_NP")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 CT", "L6_CT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 IT", "L6_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 IT Car3", "L6_IT_Car3")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "Sst Chodl", "Sst_Chodl")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "MicroPVM", "Micro_PVM")

# Subset modeling results
enrichment_results <- modeling_results[["enrichment"]]

# Load bulk data
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
bulk <- bulk[, vars]

# Overlap to get genes in both datasets
overlap <- intersect(rownames(bulk), rownames(enrichment_results))
bulk <- bulk[overlap, ]
enrichment_results <- enrichment_results[overlap, ]

# Define upregulated and downregulated DEGs for PTSD and MDD
DE_bulk_PTSD_up <- rownames(bulk[bulk$dACC_adjPVal_PTSD < 0.1 & bulk$dACC_logFC_PTSD > 0, ])
DE_bulk_PTSD_down <- rownames(bulk[bulk$dACC_adjPVal_PTSD < 0.1 & bulk$dACC_logFC_PTSD < 0, ])
DE_bulk_MDD_up <- rownames(bulk[bulk$dACC_adjPVal_MDD < 0.1 & bulk$dACC_logFC_MDD > 0, ])
DE_bulk_MDD_down <- rownames(bulk[bulk$dACC_adjPVal_MDD < 0.1 & bulk$dACC_logFC_MDD < 0, ])

nonDE_bulk_PTSD_up <- setdiff(rownames(bulk), DE_bulk_PTSD_up)
nonDE_bulk_PTSD_down <- setdiff(rownames(bulk), DE_bulk_PTSD_down)
nonDE_bulk_MDD_up <- setdiff(rownames(bulk), DE_bulk_MDD_up)
nonDE_bulk_MDD_down <- setdiff(rownames(bulk), DE_bulk_MDD_down)

# Create separate results lists for PTSD and MDD
results_PTSD_up <- list()
results_PTSD_down <- list()
results_MDD_up <- list()
results_MDD_down <- list()

# Unique cell types excluding Sst_Chodl
k <- unique(sce$cellType_azimuth)
k <- k[k != "Sst_Chodl"]

#top_n <- 100

# Loop through each cell type to perform tests
for (i in k) {
    print(i)
    # Identify upregulated and downregulated DE genes in the spatial domain

    #t_stat_threshold <- sort(enrichment_results[[paste0("t_stat_", i)]], decreasing = T)[top_n]
    #DE_clust_genes_up <- rownames(enrichment_results[enrichment_results[[paste0("t_stat_", i)]] >= t_stat_threshold, ])

    #DE_clust_genes_down <- DE_clust_genes_up

    #nonDE_clust_genes_up <- rownames(enrichment_results[enrichment_results[[paste0("t_stat_", i)]] < t_stat_threshold, ])
    #nonDE_clust_genes_down <- nonDE_clust_genes_up

    DE_clust_genes_up <- rownames(enrichment_results[
        enrichment_results[[paste0("fdr_", i)]] < 0.05 & enrichment_results[[paste0("logFC_", i)]] > 1, ])
    print(length(DE_clust_genes_up))
    DE_clust_genes_down <- DE_clust_genes_up

    nonDE_clust_genes_up <- rownames(enrichment_results[
        enrichment_results[[paste0("fdr_", i)]] >= 0.05 | enrichment_results[[paste0("logFC_", i)]] <= 1, ])
    nonDE_clust_genes_down <- nonDE_clust_genes_up

    # PTSD results
    PTSD_results <- function(DE_bulk, nonDE_bulk, DE_clust_genes, nonDE_clust_genes) {
        DE_clust_DE_bulk <- length(intersect(DE_clust_genes, DE_bulk))
        DE_clust_nonDE_bulk <- length(intersect(DE_clust_genes, nonDE_bulk))
        nonDE_clust_DE_bulk <- length(intersect(nonDE_clust_genes, DE_bulk))
        nonDE_clust_nonDE_bulk <- length(intersect(nonDE_clust_genes, nonDE_bulk))

        # Create 2x2 contingency table
        contingency_table <- matrix(c(DE_clust_DE_bulk, DE_clust_nonDE_bulk,
                                      nonDE_clust_DE_bulk, nonDE_clust_nonDE_bulk),
                                    nrow = 2, byrow = TRUE,
                                    dimnames = list("Cluster" = c("DE", "nonDE"),
                                                    "Bulk" = c("DE", "nonDE")))
        # Perform Fisher's exact test
        fisher_test <- fisher.test(contingency_table)

        list(contingency_table = contingency_table, fisher_test = fisher_test)
    }

    results_PTSD_up[[i]] <- PTSD_results(DE_bulk_PTSD_up, nonDE_bulk_PTSD_up, DE_clust_genes_up, nonDE_clust_genes_up)
    results_PTSD_down[[i]] <- PTSD_results(DE_bulk_PTSD_down, nonDE_bulk_PTSD_down, DE_clust_genes_down, nonDE_clust_genes_down)
    results_MDD_up[[i]] <- PTSD_results(DE_bulk_MDD_up, nonDE_bulk_MDD_up, DE_clust_genes_up, nonDE_clust_genes_up)
    results_MDD_down[[i]] <- PTSD_results(DE_bulk_MDD_down, nonDE_bulk_MDD_down, DE_clust_genes_down, nonDE_clust_genes_down)
}

# Initialize matrices for storing p-values
pvalues_PTSD_up <- matrix(NA, nrow = length(k), ncol = 1)
pvalues_PTSD_down <- matrix(NA, nrow = length(k), ncol = 1)
pvalues_MDD_up <- matrix(NA, nrow = length(k), ncol = 1)
pvalues_MDD_down <- matrix(NA, nrow = length(k), ncol = 1)

# Fill matrices with p-values
for (i in 1:length(k)) {
    index <- k[i]
    pvalues_PTSD_up[i, ] <- results_PTSD_up[[index]]$fisher_test$p.value
    pvalues_PTSD_down[i, ] <- results_PTSD_down[[index]]$fisher_test$p.value
    pvalues_MDD_up[i, ] <- results_MDD_up[[index]]$fisher_test$p.value
    pvalues_MDD_down[i, ] <- results_MDD_down[[index]]$fisher_test$p.value
}

# Combine matrices into one matrix for the heatmap
combined_pvalues <- cbind(
    "PTSD Up" = pvalues_PTSD_up,
    "PTSD Down" = pvalues_PTSD_down,
    "MDD Up" = pvalues_MDD_up,
    "MDD Down" = pvalues_MDD_down
)

# Assign row names for spatial domains and column names for regulation types
rownames(combined_pvalues) <- k
colnames(combined_pvalues) <- c("PTSD Up", "PTSD Down", "MDD Up", "MDD Down")

# Compute average values for reordering
row_means <- rowMeans(-log10(combined_pvalues), na.rm = TRUE)
col_means <- colMeans(-log10(combined_pvalues), na.rm = TRUE)

# Order rows and columns based on average values
ordered_rows <- order(row_means, decreasing = TRUE)
ordered_cols <- order(col_means, decreasing = TRUE)

# Reorder the heatmap matrix
combined_pvalues_ordered <- combined_pvalues[ordered_rows, ordered_cols]

col_fun <- colorRamp2(
    c(1.3, max(-log10(combined_pvalues_ordered))),
    c("white", "blue") # White for -log10(p) >= 1.3 (p >= 0.05), blue for more significant p-values
)

# Create heatmap for the combined p-values
heatmap_combined <- Heatmap(
    -log10(combined_pvalues_ordered),
    name = "-log(p)",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_side = "left",
    column_title = "snRNA-seq Psych. Enrichment",
    heatmap_legend_param = list(
        title = "-log(p)",
        title_position = "topcenter",
        title_gp = gpar(fontsize = 10),
        labels_gp = gpar(fontsize = 8)
    ),
    column_names_side = "bottom",
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 2)
)


# Display the heatmap
pdf(here("plots", "19_bulk_deconvolution", "single_nucleus_heatmap_up_down_pval_0.1_separated.pdf"), height = 4, width = 4)
draw(heatmap_combined, merge_legend = F, annotation_legend_side = "bottom")
grid.text("",
          x = unit(0.5, "npc"), y = unit(0.02, "npc"),
          just = "center", gp = gpar(fontsize = 10, col = "black"))
dev.off()
