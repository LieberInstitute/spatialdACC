setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(scuttle)
library(edgeR)
library(dplyr)

# load DEG data
load(
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_DE", "nnSVG_PRECAST_captureArea_9.Rdata")
)

# create a matrix where rows are genes, columns are clusters (domains) and elements are some statistic from DE analysis (eg t-stats)
# we will use the t-statistic for now

enrichment_results <- modeling_results[["enrichment"]]

# find repeated gene names in enrichment_results$gene
duplicated_genes <- enrichment_results[duplicated(enrichment_results$gene),]$gene

# remove one of the duplicated genes at random for each gene
for (gene in duplicated_genes) {
    enrichment_results <- enrichment_results[-sample(which(enrichment_results$gene == gene), 1),]
}

aggregated <- data.frame(
    L1 = enrichment_results$t_stat_clustL1,
    L2 = enrichment_results$t_stat_clustL2,
    L3 = enrichment_results$t_stat_clustL3,
    L5 = enrichment_results$t_stat_clustL5,
    L6a = enrichment_results$t_stat_clustL6a,
    L6b = enrichment_results$t_stat_clustL6b,
    WM = enrichment_results$t_stat_clustWM)

rownames(aggregated) <- enrichment_results$gene

write.table(aggregated, here::here("processed-data", "17-2_LDSC_enrichment", "visium_aggregated_de.tsv"), na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")
