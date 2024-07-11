setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library(here)

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
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "Micro-PVM", "Micro_PVM")


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
# [1] 17292

# subset the bulk object
bulk <- bulk[overlap,]

# subset the modeling results
enrichment_results <- enrichment_results[overlap,]

DE_bulk_genes <- rownames(bulk[bulk$dACC_adjPVal_PTSD < 0.05 | bulk$dACC_adjPVal_MDD < 0.05,])
nonDE_bulk_genes <- rownames(bulk[bulk$dACC_adjPVal_PTSD >= 0.05 & bulk$dACC_adjPVal_MDD >= 0.05,])

# create a list to store the results
results <- list()

k <- unique(sce$cellType_azimuth)
k <- k[k != "Sst_Chodl"]

for (i in k) {

    # identify DE and non-DE genes
    DE_clust_genes <- rownames(enrichment_results[enrichment_results[[paste0("fdr_", i)]] < 0.05, ])
    nonDE_clust_genes <- rownames(enrichment_results[enrichment_results[[paste0("fdr_", i)]] >= 0.05, ])

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

    # Perform the chi-square test
    chi_sq_test <- chisq.test(contingency_table)

    print(chi_sq_test)

    if (chi_sq_test$p.value < 0.05) {
        print(paste0("Cluster ", i, " is significant."))
    } else {
        print(paste0("Cluster ", i, " is not significant."))
    }

    fisher_test <- fisher.test(contingency_table)

    print(fisher_test)

    # Store the results
    results[[i]] <- list(contingency_table = contingency_table,
                         chi_sq_test = chi_sq_test,
                         fisher_test = fisher_test)

}

# Micro PVM is close to significant 0.59 with 0.05 bulk DE genes cutoff
