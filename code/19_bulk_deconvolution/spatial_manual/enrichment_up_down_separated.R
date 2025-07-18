setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library(here)
library(ComplexHeatmap)
library(circlize)

# Load DEG data
load(
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_DE", "nnSVG_PRECAST_captureArea_9.Rdata")
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

# Define upregulated and downregulated DEGs in bulk
DE_bulk_PTSD_up <- rownames(bulk[bulk$dACC_adjPVal_PTSD < 0.1 & bulk$dACC_logFC_PTSD > 0, ])
DE_bulk_PTSD_down <- rownames(bulk[bulk$dACC_adjPVal_PTSD < 0.1 & bulk$dACC_logFC_PTSD < 0, ])

DE_bulk_MDD_up <- rownames(bulk[bulk$dACC_adjPVal_MDD < 0.1 & bulk$dACC_logFC_MDD > 0, ])
DE_bulk_MDD_down <- rownames(bulk[bulk$dACC_adjPVal_MDD < 0.1 & bulk$dACC_logFC_MDD < 0, ])

# Define non-DE sets for PTSD and MDD, up and down
nonDE_bulk_PTSD_up <- setdiff(rownames(bulk), DE_bulk_PTSD_up)
nonDE_bulk_PTSD_down <- setdiff(rownames(bulk), DE_bulk_PTSD_down)
nonDE_bulk_MDD_up <- setdiff(rownames(bulk), DE_bulk_MDD_up)
nonDE_bulk_MDD_down <- setdiff(rownames(bulk), DE_bulk_MDD_down)

# Create a list to store the results for each condition and direction
results_PTSD_up <- list()
results_PTSD_down <- list()
results_MDD_up <- list()
results_MDD_down <- list()

for (i in c("L1","L2","L2_3","L4_5","L5","L6","CC","WM")) {
    # Identify upregulated and downregulated DE genes in the spatial domain
    DE_clust_genes_up <- rownames(enrichment_results[
        enrichment_results[[paste0("fdr_", i)]] < 0.05 & enrichment_results[[paste0("logFC_", i)]] > 0, ])
    DE_clust_genes_down <- rownames(enrichment_results[
        enrichment_results[[paste0("fdr_", i)]] < 0.05 & enrichment_results[[paste0("logFC_", i)]] < 0, ])

    nonDE_clust_genes_up <- rownames(enrichment_results[
        enrichment_results[[paste0("fdr_", i)]] >= 0.05 | enrichment_results[[paste0("logFC_", i)]] <= 0, ])
    nonDE_clust_genes_down <- rownames(enrichment_results[
        enrichment_results[[paste0("fdr_", i)]] >= 0.05 | enrichment_results[[paste0("logFC_", i)]] >= 0, ])

    # PTSD upregulated
    DE_clust_DE_bulk_PTSD_up <- length(intersect(DE_clust_genes_up, DE_bulk_PTSD_up))
    DE_clust_nonDE_bulk_PTSD_up <- length(intersect(DE_clust_genes_up, nonDE_bulk_PTSD_up))
    nonDE_clust_DE_bulk_PTSD_up <- length(intersect(nonDE_clust_genes_up, DE_bulk_PTSD_up))
    nonDE_clust_nonDE_bulk_PTSD_up <- length(intersect(nonDE_clust_genes_up, nonDE_bulk_PTSD_up))

    contingency_table_PTSD_up <- matrix(c(DE_clust_DE_bulk_PTSD_up, DE_clust_nonDE_bulk_PTSD_up,
                                          nonDE_clust_DE_bulk_PTSD_up, nonDE_clust_nonDE_bulk_PTSD_up),
                                        nrow = 2, byrow = TRUE,
                                        dimnames = list("Cluster" = c("DE", "nonDE"),
                                                        "Bulk" = c("DE", "nonDE")))

    fisher_test_PTSD_up <- fisher.test(contingency_table_PTSD_up)
    results_PTSD_up[[i]] <- list(fisher_test = fisher_test_PTSD_up)

    # PTSD downregulated
    DE_clust_DE_bulk_PTSD_down <- length(intersect(DE_clust_genes_down, DE_bulk_PTSD_down))
    DE_clust_nonDE_bulk_PTSD_down <- length(intersect(DE_clust_genes_down, nonDE_bulk_PTSD_down))
    nonDE_clust_DE_bulk_PTSD_down <- length(intersect(nonDE_clust_genes_down, DE_bulk_PTSD_down))
    nonDE_clust_nonDE_bulk_PTSD_down <- length(intersect(nonDE_clust_genes_down, nonDE_bulk_PTSD_down))

    contingency_table_PTSD_down <- matrix(c(DE_clust_DE_bulk_PTSD_down, DE_clust_nonDE_bulk_PTSD_down,
                                            nonDE_clust_DE_bulk_PTSD_down, nonDE_clust_nonDE_bulk_PTSD_down),
                                          nrow = 2, byrow = TRUE,
                                          dimnames = list("Cluster" = c("DE", "nonDE"),
                                                          "Bulk" = c("DE", "nonDE")))

    fisher_test_PTSD_down <- fisher.test(contingency_table_PTSD_down)
    results_PTSD_down[[i]] <- list(fisher_test = fisher_test_PTSD_down)

    # MDD upregulated
    DE_clust_DE_bulk_MDD_up <- length(intersect(DE_clust_genes_up, DE_bulk_MDD_up))
    DE_clust_nonDE_bulk_MDD_up <- length(intersect(DE_clust_genes_up, nonDE_bulk_MDD_up))
    nonDE_clust_DE_bulk_MDD_up <- length(intersect(nonDE_clust_genes_up, DE_bulk_MDD_up))
    nonDE_clust_nonDE_bulk_MDD_up <- length(intersect(nonDE_clust_genes_up, nonDE_bulk_MDD_up))

    contingency_table_MDD_up <- matrix(c(DE_clust_DE_bulk_MDD_up, DE_clust_nonDE_bulk_MDD_up,
                                         nonDE_clust_DE_bulk_MDD_up, nonDE_clust_nonDE_bulk_MDD_up),
                                       nrow = 2, byrow = TRUE,
                                       dimnames = list("Cluster" = c("DE", "nonDE"),
                                                       "Bulk" = c("DE", "nonDE")))

    fisher_test_MDD_up <- fisher.test(contingency_table_MDD_up)
    results_MDD_up[[i]] <- list(fisher_test = fisher_test_MDD_up)

    # MDD downregulated
    DE_clust_DE_bulk_MDD_down <- length(intersect(DE_clust_genes_down, DE_bulk_MDD_down))
    DE_clust_nonDE_bulk_MDD_down <- length(intersect(DE_clust_genes_down, nonDE_bulk_MDD_down))
    nonDE_clust_DE_bulk_MDD_down <- length(intersect(nonDE_clust_genes_down, DE_bulk_MDD_down))
    nonDE_clust_nonDE_bulk_MDD_down <- length(intersect(nonDE_clust_genes_down, nonDE_bulk_MDD_down))

    contingency_table_MDD_down <- matrix(c(DE_clust_DE_bulk_MDD_down, DE_clust_nonDE_bulk_MDD_down,
                                           nonDE_clust_DE_bulk_MDD_down, nonDE_clust_nonDE_bulk_MDD_down),
                                         nrow = 2, byrow = TRUE,
                                         dimnames = list("Cluster" = c("DE", "nonDE"),
                                                         "Bulk" = c("DE", "nonDE")))

    fisher_test_MDD_down <- fisher.test(contingency_table_MDD_down)
    results_MDD_down[[i]] <- list(fisher_test = fisher_test_MDD_down)
}

# Combine the p-values into matrices for each condition and regulation direction
pvalues_PTSD_up <- sapply(results_PTSD_up, function(x) x$fisher_test$p.value)
pvalues_PTSD_down <- sapply(results_PTSD_down, function(x) x$fisher_test$p.value)
pvalues_MDD_up <- sapply(results_MDD_up, function(x) x$fisher_test$p.value)
pvalues_MDD_down <- sapply(results_MDD_down, function(x) x$fisher_test$p.value)

# Combine all p-value matrices
combined_pvalues <- cbind(
    PTSD_Upregulated = pvalues_PTSD_up,
    PTSD_Downregulated = pvalues_PTSD_down,
    MDD_Upregulated = pvalues_MDD_up,
    MDD_Downregulated = pvalues_MDD_down
)

# Reorder the matrix for the heatmap based on average p-values
row_means <- rowMeans(-log10(combined_pvalues), na.rm = TRUE)
col_means <- colMeans(-log10(combined_pvalues), na.rm = TRUE)

ordered_rows <- order(row_means, decreasing = TRUE)
ordered_cols <- order(col_means, decreasing = TRUE)

combined_pvalues_ordered <- combined_pvalues[ordered_rows, ordered_cols]

# Define custom color function for -log10 transformed values
col_fun <- colorRamp2(
    c(1.3, max(-log10(combined_pvalues_ordered))),
    c("white", "blue") # White for -log10(p) >= 1.3 (p >= 0.05), blue for more significant p-values
)

# Adjust heatmap
heatmap_combined <- Heatmap(
    -log10(combined_pvalues_ordered),
    name = "-log10(Fisher's p-value)",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_side = "left",
    column_title = "Spatial Domains Manual Anno DEG Enrichment",
    heatmap_legend_param = list(
        title = "-log10(p-value)",
        title_position = "topcenter",
        title_gp = gpar(fontsize = 10),
        labels_gp = gpar(fontsize = 8)
    ),
    column_names_side = "bottom",
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 2)
)


# Display the heatmap
pdf(here("plots", "19_bulk_deconvolution", "spatial_heatmap_up_down_pval_0.1_separated_manual_anno.pdf"))
draw(heatmap_combined)
grid.text("bulk cutoff changed to pval < 0.1",
          x = unit(0.5, "npc"), y = unit(0.02, "npc"),
          just = "center", gp = gpar(fontsize = 10, col = "black"))
dev.off()

