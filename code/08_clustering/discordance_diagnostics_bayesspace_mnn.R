
setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(Seurat)
    library(SpatialExperiment)
    library(PRECAST)
    library(spatialLIBD)
    library(ggplot2)
    library(gridExtra)
    library(here)
    library(fasthplus)
    library(sessioninfo)
})

load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

#import bayesspace mnn clusters
bayesSpace_name <- paste0("bayesSpace_captureArea_", k)
spe <- cluster_import(
    spe,
    cluster_dir = here::here("processed-data", "08_clustering", "BayesSpace", "preprocess_mnn", bayesSpace_name)
)

dim(reducedDims(spe)$"pp-GLM-PCA")
# [1] 77553    50

# hpb estimate. t = pre-bootstrap sample size, D = reduced dimensions matrix, L = cluster labels, r = number of bootstrap iterations

find_t <- function(L, proportion = 0.05) {
    initial_t <- floor(length(L) * proportion)
    smallest_cluster_size <- min(table(L))
    n_labels <- length(unique(L))
    ifelse(smallest_cluster_size > (initial_t / n_labels), initial_t, smallest_cluster_size * n_labels)
}

initial_t <- find_t(L = colData(spe)[[paste0("imported_bayesSpace_captureArea_", k)]], proportion = 0.01)

cluster_prop <- table(colData(spe)[[paste0("imported_bayesSpace_captureArea_", k)]]) / ncol(spe)
#         1          2          3          4          5
#0.02892216 0.76436759 0.13202584 0.05846324 0.01622116

bad_clusters <- which(cluster_prop < 0.01 / k)
if (length(bad_clusters) > 0) {
    message("For k: ", k, " we are dropping small clusters: ", paste(names(bad_clusters), collapse = ", "))
    spe <- spe[, !colData(spe)[[paste0("imported_bayesSpace_captureArea_", k)]] %in% as.integer(names(bad_clusters))]
    updated_t <- find_t(colData(spe)[[paste0("imported_bayesSpace_captureArea_", k)]], 0.01)
    message("initial t: ", initial_t, "; updated t: ", updated_t)
} else{
    updated_t <- initial_t
}

set.seed(7)
fasthplus <- hpb(D = reducedDims(spe)$"pp-GLM-PCA", L = colData(spe)[[paste0("imported_bayesSpace_captureArea_", k)]], t = updated_t, r = 30)
results <- data.frame(k = k, fasthplus = fasthplus)
write.table(results, file = here::here("processed-data", "08_clustering", "cluster_diagnostics", "fasthplus", "fasthplus_results_bayesSpace_mnn.csv"),
append = TRUE,
row.names = FALSE,
col.names =!file.exists(here("processed-data", "08_clustering", "cluster_diagnostics", "fasthplus", "fasthplus_results_bayesspace_mnn.csv")))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
