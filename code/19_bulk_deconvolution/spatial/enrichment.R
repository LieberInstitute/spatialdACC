setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library(here)

# load DEG data
load(
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_DE", "nnSVG_PRECAST_captureArea_9.Rdata")
)

# subset modeling results
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

# remove .xx from gene names
rownames(bulk) <- base_gene_ids

# list columns of interest from bulk data
vars <- c("dACC_logFC_PTSD", "dACC_t_PTSD", "dACC_adjPVal_PTSD",
          "dACC_logFC_MDD", "dACC_t_MDD", "dACC_adjPVal_MDD")

# create smaller df with these vars
bulk <- bulk[, vars]

# overlap to get genes in both datasets
overlap <- intersect(rownames(bulk), rownames(enrichment_results))
length(overlap)
# [1] 12203

# subset the bulk object
bulk <- bulk[overlap,]

# subset the modeling results
enrichment_results <- enrichment_results[overlap,]

DE_bulk_genes <- rownames(bulk[bulk$dACC_adjPVal_PTSD < 0.05 | bulk$dACC_adjPVal_MDD < 0.05,])
nonDE_bulk_genes <- rownames(bulk[bulk$dACC_adjPVal_PTSD >= 0.05 & bulk$dACC_adjPVal_MDD >= 0.05,])

# within PTSD, count how many DE genes
DE_bulk_PTSD_genes <- rownames(bulk[bulk$dACC_adjPVal_PTSD < 0.05,])
# how many of these 12 genes have logFC > 0
sum(bulk[DE_bulk_PTSD_genes, "dACC_logFC_PTSD"] <0)

# within MDD, count how many DE genes
DE_bulk_MDD_genes <- rownames(bulk[bulk$dACC_adjPVal_MDD < 0.05,])
# how many of these 52 genes have logFC > 0
sum(bulk[DE_bulk_MDD_genes, "dACC_logFC_MDD"] <0)

# which gene is overlapping between PTSD and MDD
intersect(DE_bulk_PTSD_genes, DE_bulk_MDD_genes)
# what are the fold changes for this gene ENSG00000249436
bulk["ENSG00000249436", c("dACC_logFC_PTSD", "dACC_logFC_MDD")]

# create a list to store the results
results <- list()

for (i in 1:9) {

    # identify DE and non-DE genes
    DE_clust_genes <- rownames(enrichment_results[enrichment_results[[paste0("fdr_clust", i)]] < 0.05, ])
    nonDE_clust_genes <- rownames(enrichment_results[enrichment_results[[paste0("fdr_clust", i)]] >= 0.05, ])

    # Count overlaps
    # DE in cluster and DE in bulk
    DE_clust_DE_bulk <- length(intersect(DE_clust_genes, DE_bulk_genes))

    # DE in cluster and non-DE in bulk
    DE_clust_nonDE_bulk <- length(intersect(DE_clust_genes, nonDE_bulk_genes))

    # Non-DE in cluster and DE in bulk
    nonDE_clust_DE_bulk <- length(intersect(nonDE_clust_genes, DE_bulk_genes))

    # Non-DE in cluster and non-DE in bulk
    nonDE_clust_nonDE_bulk <- length(intersect(nonDE_clust_genes, nonDE_bulk_genes))

    # Create the 2x2 table
    contingency_table <- matrix(c(DE_clust_DE_bulk, DE_clust_nonDE_bulk,
                                  nonDE_clust_DE_bulk, nonDE_clust_nonDE_bulk),
                                nrow = 2, byrow = TRUE,
                                dimnames = list("Cluster" = c("DE", "nonDE"),
                                                "Bulk" = c("DE", "nonDE")))

    print(contingency_table)

    # Perform the chi-square test
    chi_sq_test <- chisq.test(contingency_table)

    print(chi_sq_test)

    fisher_test <- fisher.test(contingency_table)

    print(fisher_test)

    # Store the results
    results[[i]] <- list(contingency_table = contingency_table,
                         chi_sq_test = chi_sq_test,
                         fisher_test = fisher_test)

}

# clusters 3 (WM1) and 9 (L1) are significant using 0.05 for bulk DE cutoff

#print all chi  sq tests
for (i in 1:9) {
    print(results[[i]]$fisher_test)
}
