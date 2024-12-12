setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(tidyverse)

loads <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_patterns_subset.RDS"))

no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

ldsc.score <- as.matrix(read.csv(here::here("processed-data", "17_LDSC", "NMF_score.csv"),
                                 row.names = 1))

# make the last underscore - in colnames(ldsc.score) to .
colnames(loads) <- gsub("-", ".", colnames(loads))

colnames(ldsc.score)[25] <- "SST_Chodl.NMF51"

#### new nmf_score_top-rank.csv
identical(rownames(loads), rownames(ldsc.score))
identical(colnames(loads), colnames(ldsc.score))
topN.mat = matrix(0, nrow=nrow(ldsc.score), ncol=ncol(ldsc.score))


rownames(topN.mat) <- rownames(ldsc.score)
colnames(topN.mat) <- colnames(ldsc.score)

for(i in colnames(ldsc.score)) {
    tmp = loads[,i]
    rank1 = (1+length(tmp))-rank(tmp)
    tmp.bin = ifelse(rank1<=928, 1, 0)
    stopifnot(identical(names(tmp.bin), rownames(topN.mat)))
    topN.mat[,i] = tmp.bin
}

saveRDS(topN.mat, here("processed-data", "17_LDSC", "top928_intermediate_mat.rds"))
