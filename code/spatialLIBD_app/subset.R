library("spatialLIBD")
library("lobstr")
library("here")
library("sessioninfo")


load("modeling-nnSVG_PRECAST_captureArea_9.Rdata", verbose = TRUE)
load("nnSVG_PRECAST_captureArea_9.Rdata", verbose = TRUE)

# Rscript to subset for sig_genes

sig_genes <- sig_genes_extract_all(
    n = nrow(spe_pseudo),
    modeling_results = modeling_results,
    sce_layer = spe_pseudo
)

lobstr::obj_size(sig_genes)
# 205.76 MB

dim(sig_genes)
# [1] [1] 678800     13

## Drop parts we don't need to reduce the memory
sig_genes$in_rows <- NULL
sig_genes$in_rows_top20 <- NULL
lobstr::obj_size(sig_genes)
# 59.13 MB

save(sig_genes,file = "sig_genes_subset.Rdata")
