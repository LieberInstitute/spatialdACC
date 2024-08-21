setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(RcppML)
library(here)

x <- readRDS(file = here("processed-data", "18_PsychENCODE_NMF", "nmf_results.rds"))
sce <- readRDS(file = here("processed-data", "18_PsychENCODE_NMF", "pseudobulk_combined.rds"))

# function for getting top n genes for each pattern
top_genes <- function(W, n=10){
    top_genes <- apply(W, 2, function(x) names(sort(x, decreasing=TRUE)[1:n]))
    return(top_genes)
}

# get top 10 genes
top10 <- top_genes(x$w, 10)
write.csv(top10, file = here("processed-data", "18_PsychENCODE_NMF", "top10_genes.csv"))
