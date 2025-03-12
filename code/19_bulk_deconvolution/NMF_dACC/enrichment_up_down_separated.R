setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library(here)
library(ComplexHeatmap)
library(circlize)
library(SingleCellExperiment)
library(grid)
library(SpatialExperiment)
library(RcppML)

x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))
load(file = here("processed-data", "snRNA-seq", "03_batch_correction", "sce_harmony.Rdata"))

# function for getting top n genes for each pattern
top_genes <- function(W, n=10){
    top_genes <- apply(W, 2, function(x) names(sort(x, decreasing=TRUE)[1:n]))
    return(top_genes)
}

# get top 500 genes
top500 <- top_genes(x$w, 500)

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

# need to get all bg genes from NMF as well
bg_genes_id <- rowData(sce)$gene_id

# Overlap to get genes in both datasets
overlap <- intersect(rownames(bulk), bg_genes_id)
bulk <- bulk[overlap, ]
bg_genes_id <- intersect(bg_genes_id,overlap)

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

gene_list <- list()
gene_list[["nmf38"]] <- setdiff(top500[,"nmf38"],top500[,"nmf61"])
gene_list[["nmf61"]] <- setdiff(top500[,"nmf61"],top500[,"nmf38"])

k <- c("nmf38","nmf61")

for (i in k) {
    print(i)

    top500_i <- gene_list[[i]]

    # match gene names to gene ids from sce
    top500_i_id <- rowData(sce)$gene_id[match(top500_i, rowData(sce)$gene_name)]

    # make sure top genes are in the intersection
    top500_i_id <- intersect(top500_i_id,overlap)

    #background genes
    bg_genes_i_id <- setdiff(bg_genes_id,top500_i_id)

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

    results_PTSD_up[[i]] <- PTSD_results(DE_bulk_PTSD_up, nonDE_bulk_PTSD_up, top500_i_id, bg_genes_i_id)
    results_PTSD_down[[i]] <- PTSD_results(DE_bulk_PTSD_down, nonDE_bulk_PTSD_down, top500_i_id, bg_genes_i_id)
    results_MDD_up[[i]] <- PTSD_results(DE_bulk_MDD_up, nonDE_bulk_MDD_up, top500_i_id, bg_genes_i_id)
    results_MDD_down[[i]] <- PTSD_results(DE_bulk_MDD_down, nonDE_bulk_MDD_down, top500_i_id, bg_genes_i_id)
}
