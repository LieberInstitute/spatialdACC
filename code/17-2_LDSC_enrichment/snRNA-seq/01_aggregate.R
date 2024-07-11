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

enrichment_results <- modeling_results[["enrichment"]]

which(enrichment_results$gene == "HSPA14")
#[1]  10602 15447
enrichment_results <- enrichment_results[-10602,]

# find repeated gene names in enrichment_results$gene
duplicated_genes <- enrichment_results[duplicated(enrichment_results$gene),]$gene

# remove one of the duplicated genes at random for each gene
for (gene in duplicated_genes) {
    enrichment_results <- enrichment_results[-sample(which(enrichment_results$gene == gene), 1),]
}

k <- unique(sce$cellType_azimuth)
k <- k[k != "Sst_Chodl"]

aggregated <- enrichment_results[,c(1:length(k))]

#remove "t_stat_" prefix from column names
colnames(aggregated) <- gsub("t_stat_", "", colnames(aggregated))

rownames(aggregated) <- enrichment_results$gene

write.table(aggregated, here::here("processed-data", "17-2_LDSC_enrichment", "snRNA-seq_aggregated_de.tsv"), na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")
