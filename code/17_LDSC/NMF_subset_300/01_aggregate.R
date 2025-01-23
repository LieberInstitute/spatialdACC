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

# remove repeated patterns
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "misc.NMF15")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "misc.NMF59")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "misc.NMF61")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "VLMC.NMF17")]

# choose only one pattern for each cell type/class
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Oligo.NMF26")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Oligo.NMF23")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Oligo.NMF13")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Oligo.NMF43")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Oligo.NMF40")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Oligo.NMF36")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Oligo.NMF28")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Oligo.NMF9")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Oligo.NMF33")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Oligo.NMF39")]

ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "L2_3_IT.NMF11")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Pvalb.NMF63")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "SST.NMF56")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "LAMP5.NMF37")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Sncg.NMF58")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Vip.NMF47")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Endo.NMF49")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Astro.NMF14")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Astro.NMF53")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "Astro.NMF65")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "OPC.NMF24")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "microPVM.NMF54")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "microPVM.NMF57")]
ldsc.score <- ldsc.score[, !(colnames(ldsc.score) == "misc.NMF64")]

#### new nmf_score_top-rank.csv
topN.mat = matrix(0, nrow=nrow(ldsc.score), ncol=ncol(ldsc.score))

rownames(topN.mat) <- rownames(ldsc.score)
colnames(topN.mat) <- colnames(ldsc.score)

for(i in colnames(ldsc.score)) {
    tmp = loads[,i]
    rank1 = (1+length(tmp))-rank(tmp)
    tmp.bin = ifelse(rank1<=300, 1, 0)
    stopifnot(identical(names(tmp.bin), rownames(topN.mat)))
    topN.mat[,i] = tmp.bin
}

saveRDS(topN.mat, here("processed-data", "17_LDSC", "top300_subset_intermediate_mat.rds"))
