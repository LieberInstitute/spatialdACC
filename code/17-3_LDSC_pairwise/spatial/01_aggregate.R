setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(scuttle)
library(edgeR)
library(dplyr)

load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9.Rdata"))

spe$PRECAST_cluster <- unfactor(spe$PRECAST_cluster)
spe$PRECAST_cluster[spe$PRECAST_cluster == 3] <- "WM1"
spe$PRECAST_cluster[spe$PRECAST_cluster == 8] <- "WM2"
spe$PRECAST_cluster[spe$PRECAST_cluster == 7] <- "WM-CC"
spe$PRECAST_cluster[spe$PRECAST_cluster == 5] <- "L6b"
spe$PRECAST_cluster[spe$PRECAST_cluster == 6] <- "L6a"
spe$PRECAST_cluster[spe$PRECAST_cluster == 4] <- "L5"
spe$PRECAST_cluster[spe$PRECAST_cluster == 2] <- "L3"
spe$PRECAST_cluster[spe$PRECAST_cluster == 1] <- "L2"
spe$PRECAST_cluster[spe$PRECAST_cluster == 9] <- "L1"


# load DEG data
load(
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_DE", "nnSVG_PRECAST_captureArea_9.Rdata")
)

# create a matrix where rows are genes, columns are clusters (domains) and elements are some statistic from DE analysis (eg t-stats)
# we will use the t-statistic for now

pairwise_results <- modeling_results[["pairwise"]]

# find repeated gene names in pairwise_results$gene
duplicated_genes <- pairwise_results[duplicated(pairwise_results$gene),]$gene

# remove one of the duplicated genes at random for each gene
for (gene in duplicated_genes) {
    pairwise_results <- pairwise_results[-sample(which(pairwise_results$gene == gene), 1),]
}

gene <- pairwise_results$gene

# keep columns 1-36
pairwise_results <- pairwise_results[,1:36]
colnames(pairwise_results) <- gsub("t_stat_", "", colnames(pairwise_results))

# initialize a matrix to store aggregated results (rows = genes, columns = clusters 1-9)
aggregated <- matrix(NA, nrow = nrow(pairwise_results), ncol = 9)
# make col names "clust1" to "clust9"
colnames(aggregated) <- paste0("clust", 1:9)

# loop through each COLUMN in pairwise_results
cols <- colnames(pairwise_results)

for (c in cols) {
    # get the cluster number
    first_clust <- strsplit(c, "-")[[1]][1]
    second_clust <- strsplit(c, "-")[[1]][2]

    # loop through each gene for this column
    for (i in 1:nrow(pairwise_results)) {
        # Get the t-statistic for this gene and cluster
        t_stat <- pairwise_results[i, c]

        # pf the t-statistic is positive, add to the first cluster, if negative, add to the second cluster as a positive
        if (t_stat > 0) {
            if (is.na(aggregated[i, first_clust]) || t_stat > aggregated[i, first_clust]) {
                aggregated[i, first_clust] <- t_stat
            }
        } else {
            if (is.na(aggregated[i, second_clust]) || -t_stat > aggregated[i, second_clust]) {
                aggregated[i, second_clust] <- -t_stat
            }
        }
    }
}

aggregated <- as.data.frame(aggregated)
rownames(aggregated) <- gene

write.table(aggregated, here::here("processed-data", "17-3_LDSC_pairwise", "visium_aggregated_de.tsv"), na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")
