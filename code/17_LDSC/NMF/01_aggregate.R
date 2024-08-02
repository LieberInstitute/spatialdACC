setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(scuttle)
library(edgeR)

x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))

patterns <- x@w

# remove NMF factors related to batch or technical vars
# 66, 1, 23, 8, 21
patterns <- patterns[, -c(66, 1, 23, 8, 21)]

# remove 16,134 genes that have all 0s
all_zeros <- which(rowSums(patterns)==0)
patterns <- patterns[-all_zeros,]

# calculate z-scores for each gene
scaled <- t(scale(t(patterns)))

write.table(scaled, here::here("processed-data", "17_LDSC", "NMF_aggregated_cpm.tsv"), na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")

