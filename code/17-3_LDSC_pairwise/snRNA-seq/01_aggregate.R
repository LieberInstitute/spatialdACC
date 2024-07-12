setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(scuttle)
library(edgeR)
library(dplyr)

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))
sce <- logNormCounts(sce)

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

# load DEG data
load(
    file = here("processed-data", "snRNA-seq", "05_azimuth", "azimuth_DE_results.Rdata")
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

k <- unique(sce$cellType_azimuth)
k <- k[k != "Sst_Chodl"]

# keep columns 1-36
pairwise_results <- pairwise_results[,1:171]
colnames(pairwise_results) <- gsub("t_stat_", "", colnames(pairwise_results))

# initialize a matrix to store aggregated results (rows = genes, columns = clusters)
aggregated <- matrix(NA, nrow = nrow(pairwise_results), ncol = length(k))
colnames(aggregated) <- k

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

write.table(aggregated, here::here("processed-data", "17-3_LDSC_pairwise", "snRNA-seq_aggregated_de.tsv"), na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")
