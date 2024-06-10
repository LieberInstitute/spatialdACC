setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(RcppML)
library(here)
library(spatialLIBD)

sce_path_zip <- fetch_data("spatialDLPFC_snRNAseq")
sce_path <- unzip(sce_path_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
    file.path(tempdir(), "sce_DLPFC_annotated")
)
x <- readRDS(file = here("processed-data", "15_cross_region_snRNA-seq", "DLPFC_nmf_results.RDS"))

# function for getting top n genes for each pattern
top_genes <- function(W, n=10){
    top_genes <- apply(W, 2, function(x) names(sort(x, decreasing=TRUE)[1:n]))
    return(top_genes)
}

# get top 10 genes
top10 <- top_genes(x$w, 10)
write.csv(top10, file = here("processed-data", "15_cross_region_snRNA-seq", "top10_genes.csv"))

# get top 20 genes
top20 <- top_genes(x$w, 20)
write.csv(top20, file = here("processed-data", "15_cross_region_snRNA-seq", "top20_genes.csv"))
