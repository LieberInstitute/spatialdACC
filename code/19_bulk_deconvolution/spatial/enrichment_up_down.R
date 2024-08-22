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
DE_bulk_PTSD_up <- rownames(bulk[bulk$dACC_adjPVal_PTSD < 0.05 & bulk$dACC_logFC_PTSD > 0, ])
DE_bulk_PTSD_down <- rownames(bulk[bulk$dACC_adjPVal_PTSD < 0.05 & bulk$dACC_logFC_PTSD < 0, ])

DE_bulk_MDD_up <- rownames(bulk[bulk$dACC_adjPVal_MDD < 0.05 & bulk$dACC_logFC_MDD > 0, ])
DE_bulk_MDD_down <- rownames(bulk[bulk$dACC_adjPVal_MDD < 0.05 & bulk$dACC_logFC_MDD < 0, ])

# Combine PTSD and MDD DEGs for upregulated and downregulated separately
DE_bulk_up <- unique(c(DE_bulk_PTSD_up, DE_bulk_MDD_up))
DE_bulk_down <- unique(c(DE_bulk_PTSD_down, DE_bulk_MDD_down))

nonDE_bulk_up <- setdiff(rownames(bulk), DE_bulk_up)
nonDE_bulk_down <- setdiff(rownames(bulk), DE_bulk_down)

# Create a list to store the results
results_up <- list()
results_down <- list()

for (i in 1:9) {

    # Identify upregulated and downregulated DE genes in the spatial domain
    DE_clust_genes_up <- rownames(enrichment_results[
        enrichment_results[[paste0("fdr_clust", i)]] < 0.05 & enrichment_results[[paste0("logFC_clust", i)]] > 0, ])
    DE_clust_genes_down <- rownames(enrichment_results[
        enrichment_results[[paste0("fdr_clust", i)]] < 0.05 & enrichment_results[[paste0("logFC_clust", i)]] < 0, ])

    nonDE_clust_genes_up <- rownames(enrichment_results[
        enrichment_results[[paste0("fdr_clust", i)]] >= 0.05 | enrichment_results[[paste0("logFC_clust", i)]] <= 0, ])
    nonDE_clust_genes_down <- rownames(enrichment_results[
        enrichment_results[[paste0("fdr_clust", i)]] >= 0.05 | enrichment_results[[paste0("logFC_clust", i)]] >= 0, ])

    # Count overlaps for upregulated genes
    DE_clust_DE_bulk_up <- length(intersect(DE_clust_genes_up, DE_bulk_up))
    DE_clust_nonDE_bulk_up <- length(intersect(DE_clust_genes_up, nonDE_bulk_up))
    nonDE_clust_DE_bulk_up <- length(intersect(nonDE_clust_genes_up, DE_bulk_up))
    nonDE_clust_nonDE_bulk_up <- length(intersect(nonDE_clust_genes_up, nonDE_bulk_up))

    # Create the 2x2 table for upregulated genes
    contingency_table_up <- matrix(c(DE_clust_DE_bulk_up, DE_clust_nonDE_bulk_up,
                                     nonDE_clust_DE_bulk_up, nonDE_clust_nonDE_bulk_up),
                                   nrow = 2, byrow = TRUE,
                                   dimnames = list("Cluster" = c("DE", "nonDE"),
                                                   "Bulk" = c("DE", "nonDE")))

    # Perform the chi-square test and Fisher's exact test for upregulated genes
    chi_sq_test_up <- chisq.test(contingency_table_up)
    fisher_test_up <- fisher.test(contingency_table_up)

    # Store the results for upregulated genes
    results_up[[i]] <- list(contingency_table = contingency_table_up,
                            chi_sq_test = chi_sq_test_up,
                            fisher_test = fisher_test_up)

    # Count overlaps for downregulated genes
    DE_clust_DE_bulk_down <- length(intersect(DE_clust_genes_down, DE_bulk_down))
    DE_clust_nonDE_bulk_down <- length(intersect(DE_clust_genes_down, nonDE_bulk_down))
    nonDE_clust_DE_bulk_down <- length(intersect(nonDE_clust_genes_down, DE_bulk_down))
    nonDE_clust_nonDE_bulk_down <- length(intersect(nonDE_clust_genes_down, nonDE_bulk_down))

    # Create the 2x2 table for downregulated genes
    contingency_table_down <- matrix(c(DE_clust_DE_bulk_down, DE_clust_nonDE_bulk_down,
                                       nonDE_clust_DE_bulk_down, nonDE_clust_nonDE_bulk_down),
                                     nrow = 2, byrow = TRUE,
                                     dimnames = list("Cluster" = c("DE", "nonDE"),
                                                     "Bulk" = c("DE", "nonDE")))

    # Perform the chi-square test and Fisher's exact test for downregulated genes
    chi_sq_test_down <- chisq.test(contingency_table_down)
    fisher_test_down <- fisher.test(contingency_table_down)

    # Store the results for downregulated genes
    results_down[[i]] <- list(contingency_table = contingency_table_down,
                              chi_sq_test = chi_sq_test_down,
                              fisher_test = fisher_test_down)
}

# Print results for upregulated genes
for (i in 1:9) {
    print(paste("Cluster", i, "- Upregulated Genes:"))
    print(results_up[[i]]$fisher_test)
}

# Print results for downregulated genes
for (i in 1:9) {
    print(paste("Cluster", i, "- Downregulated Genes:"))
    print(results_down[[i]]$fisher_test)
}

# Initialize matrices for upregulated and downregulated p-values
pvalues_up <- matrix(NA, nrow = 9, ncol = 1)
pvalues_down <- matrix(NA, nrow = 9, ncol = 1)

# Fill matrices with p-values from Fisher's test
for (i in 1:9) {
    pvalues_up[i, ] <- results_up[[i]]$fisher_test$p.value
    pvalues_down[i, ] <- results_down[[i]]$fisher_test$p.value
}

# Combine matrices into one matrix for the heatmap
combined_pvalues <- cbind(
    Upregulated = pvalues_up,
    Downregulated = pvalues_down
)

# Assign row names for spatial domains and column names for regulation types
rownames(combined_pvalues) <- paste0("Cluster_", 1:9)
colnames(combined_pvalues) <- c("Upregulated", "Downregulated")

# rename rownames according to
# 3- WM1
# 8- WM2
# 7- WM-CC
# 5- L6b
# 6- L6a
# 4- L5
# 2- L3
# 1- L2
# 9- L1

rownames(combined_pvalues) <- c("L2", "L3", "WM1", "L5", "L6b", "L6a", "WM-CC", "WM2", "L1")

# Compute average values for reordering
row_means <- rowMeans(-log10(combined_pvalues), na.rm = TRUE)
col_means <- colMeans(-log10(combined_pvalues), na.rm = TRUE)

# Order rows and columns based on average values
ordered_rows <- order(row_means, decreasing = TRUE)
ordered_cols <- order(col_means, decreasing = TRUE)

# Reorder the heatmap matrix
combined_pvalues_ordered <- combined_pvalues[ordered_rows, ordered_cols]


# Create heatmap for the combined p-values
heatmap_combined <- Heatmap(
    -log10(combined_pvalues_ordered),
    name = "-log10(Fisher's p-value)",
    col = colorRampPalette(c("white", "blue", "red"))(50),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_side = "left",
    column_title = "Spatial Domains DEG Enrichment",
    heatmap_legend_param = list(
        title = "-log10(p-value)",
        title_position = "topcenter", # Corrected positioning
        title_gp = gpar(fontsize = 10),
        labels_gp = gpar(fontsize = 8)
    ),
    column_names_side = "bottom",
    border = TRUE,  # Add a border around the heatmap
    border_gp = gpar(col = "black", lwd = 2)
)
# Display the heatmap
pdf(here("plots", "19_bulk_deconvolution", "spatial_heatmap_up_down.pdf"))
draw(heatmap_combined, merge_legend = F, annotation_legend_side = "bottom")
dev.off()
