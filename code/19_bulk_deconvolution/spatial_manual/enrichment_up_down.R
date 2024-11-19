setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library(here)
library(ComplexHeatmap)
library(circlize)

# Load DEG data
modeling_results <- readRDS(
    file = here("processed-data", "20_WM_comparisons", "modeling_results_WM_vs_CC.RDS")
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
bulk <- bulk[overlap, ]
enrichment_results <- enrichment_results[overlap, ]

generate_spatial_heatmap <- function(adj_pval_threshold = 0.1) {

    # Define upregulated and downregulated DEGs in bulk
    DE_bulk_PTSD_up <- rownames(bulk[bulk$dACC_adjPVal_PTSD < adj_pval_threshold & bulk$dACC_logFC_PTSD > 0, ])
    DE_bulk_PTSD_down <- rownames(bulk[bulk$dACC_adjPVal_PTSD < adj_pval_threshold & bulk$dACC_logFC_PTSD < 0, ])

    DE_bulk_MDD_up <- rownames(bulk[bulk$dACC_adjPVal_MDD < adj_pval_threshold & bulk$dACC_logFC_MDD > 0, ])
    DE_bulk_MDD_down <- rownames(bulk[bulk$dACC_adjPVal_MDD < adj_pval_threshold & bulk$dACC_logFC_MDD < 0, ])

    DE_bulk_up <- unique(c(DE_bulk_PTSD_up, DE_bulk_MDD_up))
    DE_bulk_down <- unique(c(DE_bulk_PTSD_down, DE_bulk_MDD_down))

    nonDE_bulk_up <- setdiff(rownames(bulk), DE_bulk_up)
    nonDE_bulk_down <- setdiff(rownames(bulk), DE_bulk_down)

    results_up <- list()
    results_down <- list()

    for (i in c("L1","L2","L2_3","L4_5","L5","L6","CC","WM")) {

        #check if i is in the enrichment results
        if (!paste0("fdr_", i) %in% colnames(enrichment_results)) {
            print(paste0("No enrichment results for ", i))
        }

        DE_clust_genes_up <- rownames(enrichment_results[
            enrichment_results[[paste0("fdr_", i)]] < 0.05 & enrichment_results[[paste0("logFC_", i)]] > 0, ])
        DE_clust_genes_down <- rownames(enrichment_results[
            enrichment_results[[paste0("fdr_", i)]] < 0.05 & enrichment_results[[paste0("logFC_", i)]] < 0, ])

        nonDE_clust_genes_up <- rownames(enrichment_results[
            enrichment_results[[paste0("fdr_", i)]] >= 0.05 | enrichment_results[[paste0("logFC_", i)]] <= 0, ])
        nonDE_clust_genes_down <- rownames(enrichment_results[
            enrichment_results[[paste0("fdr_", i)]] >= 0.05 | enrichment_results[[paste0("logFC_", i)]] >= 0, ])

        DE_clust_DE_bulk_up <- length(intersect(DE_clust_genes_up, DE_bulk_up))
        DE_clust_nonDE_bulk_up <- length(intersect(DE_clust_genes_up, nonDE_bulk_up))
        nonDE_clust_DE_bulk_up <- length(intersect(nonDE_clust_genes_up, DE_bulk_up))
        nonDE_clust_nonDE_bulk_up <- length(intersect(nonDE_clust_genes_up, nonDE_bulk_up))

        contingency_table_up <- matrix(c(DE_clust_DE_bulk_up, DE_clust_nonDE_bulk_up,
                                         nonDE_clust_DE_bulk_up, nonDE_clust_nonDE_bulk_up),
                                       nrow = 2, byrow = TRUE,
                                       dimnames = list("Cluster" = c("DE", "nonDE"),
                                                       "Bulk" = c("DE", "nonDE")))

        print(contingency_table_up)

        fisher_test_up <- fisher.test(contingency_table_up)
        results_up[[i]] <- list(fisher_test = fisher_test_up)

        DE_clust_DE_bulk_down <- length(intersect(DE_clust_genes_down, DE_bulk_down))
        DE_clust_nonDE_bulk_down <- length(intersect(DE_clust_genes_down, nonDE_bulk_down))
        nonDE_clust_DE_bulk_down <- length(intersect(nonDE_clust_genes_down, DE_bulk_down))
        nonDE_clust_nonDE_bulk_down <- length(intersect(nonDE_clust_genes_down, nonDE_bulk_down))

        contingency_table_down <- matrix(c(DE_clust_DE_bulk_down, DE_clust_nonDE_bulk_down,
                                           nonDE_clust_DE_bulk_down, nonDE_clust_nonDE_bulk_down),
                                         nrow = 2, byrow = TRUE,
                                         dimnames = list("Cluster" = c("DE", "nonDE"),
                                                         "Bulk" = c("DE", "nonDE")))

        fisher_test_down <- fisher.test(contingency_table_down)
        results_down[[i]] <- list(fisher_test = fisher_test_down)
    }

    pvalues_up <- sapply(results_up, function(x) x$fisher_test$p.value)
    pvalues_down <- sapply(results_down, function(x) x$fisher_test$p.value)

    combined_pvalues <- cbind(
        Upregulated = pvalues_up,
        Downregulated = pvalues_down
    )

    row_means <- rowMeans(-log10(combined_pvalues), na.rm = TRUE)
    col_means <- colMeans(-log10(combined_pvalues), na.rm = TRUE)

    ordered_rows <- order(row_means, decreasing = TRUE)
    ordered_cols <- order(col_means, decreasing = TRUE)

    combined_pvalues_ordered <- combined_pvalues[ordered_rows, ordered_cols]

    heatmap_combined <- Heatmap(
        -log10(combined_pvalues_ordered),
        name = "-log10(Fisher's p-value)",
        col = colorRampPalette(c("white", "blue", "red"))(50),
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

    print(combined_pvalues_ordered)

    # Return the heatmap object
    return(heatmap_combined)
}

# Display the heatmap
pdf(here("plots", "19_bulk_deconvolution", "spatial_heatmap_up_down_pval_0.1_manual_anno.pdf"))
generate_spatial_heatmap()
grid.text("bulk cutoff changed to pval < 0.1",
          x = unit(0.5, "npc"), y = unit(0.02, "npc"),
          just = "center", gp = gpar(fontsize = 10, col = "black"))
dev.off()
