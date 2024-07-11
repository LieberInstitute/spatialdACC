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

enrichment_results <- modeling_results[["enrichment"]]

# find repeated gene names in enrichment_results$gene
duplicated_genes <- enrichment_results[duplicated(enrichment_results$gene),]$gene

# remove one of the duplicated genes at random for each gene
for (gene in duplicated_genes) {
    enrichment_results <- enrichment_results[-sample(which(enrichment_results$gene == gene), 1),]
}

aggregated <- data.frame(
    L1 = enrichment_results$t_stat_clust9,
    L2 = enrichment_results$t_stat_clust1,
    L3 = enrichment_results$t_stat_clust2,
    L5 = enrichment_results$t_stat_clust4,
    L6a = enrichment_results$t_stat_clust6,
    L6b = enrichment_results$t_stat_clust5,
    WM1 = enrichment_results$t_stat_clust3,
    WM2 = enrichment_results$t_stat_clust8,
    WM_CC = enrichment_results$t_stat_clust7
)

rownames(aggregated) <- enrichment_results$gene

write.table(aggregated, here::here("processed-data", "17-2_LDSC_enrichment", "visium_aggregated_de.tsv"), na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")
