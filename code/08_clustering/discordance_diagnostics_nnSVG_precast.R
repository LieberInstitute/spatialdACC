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

k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

read_barcoded_csv <- function(x) {
    df <- read.csv(x)
    colnames(df) <- tolower(colnames(df))

    if (colnames(df)[2] == "cluster") {
        colnames(df)[2] <-
            gsub("gene_expression_", "", basename(dirname(x)))
    }
    return(df)
}

nnSVG_PRECAST_import <- function(spe, cluster_dir = file.path(tempdir(), "exported_clusters")) {
    clustering_files <-
        list.files(
            here::here("processed-data", "08_clustering", "PRECAST", nnSVG_precast_name),
            pattern = "clusters.csv",
            all.files = TRUE,
            full.names = TRUE,
            recursive = TRUE
        )
    clusters_list <- lapply(clustering_files, read_barcoded_csv)
    clusters <- Reduce(function(...) merge(..., by = "key", all = TRUE), clusters_list)
    cluster_cols <- which(colnames(clusters) != "key")
    colnames(clusters)[cluster_cols] <- paste0("", colnames(clusters)[cluster_cols])

    colData(spe)[,nnSVG_precast_name] <- clusters[,2]
    return(spe)
}

#import nnSVG precast clusters
nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", k)
spe <- nnSVG_PRECAST_import(
    spe,
    cluster_dir = here::here("processed-data", "08_clustering", "PRECAST", nnSVG_precast_name)
)

dim(reducedDims(spe)$"pp-GLM-PCA")
# [1] 77553    50

#find 13 spots that PRECAST filtered out
non_na_indices <- !is.na(colData(spe)[nnSVG_precast_name])
print(sum(non_na_indices))

clustering_no_na <- colData(spe)[,nnSVG_precast_name][non_na_indices]
reduced_dim_no_na <- reducedDim(spe, "pp-GLM-PCA")[non_na_indices,]

# hpb estimate. t = pre-bootstrap sample size, D = reduced dimensions matrix, L = cluster labels, r = number of bootstrap iterations

find_t <- function(L, proportion = 0.05) {
    initial_t <- floor(length(L) * proportion)
    smallest_cluster_size <- min(table(L))
    n_labels <- length(unique(L))
    ifelse(smallest_cluster_size > (initial_t / n_labels), initial_t, smallest_cluster_size * n_labels)
}

initial_t <- find_t(L = clustering_no_na, proportion = 0.01)

cluster_prop <- table(clustering_no_na) / ncol(spe)
#         1          2          3          4          5
#0.29094942 0.09810065 0.27101466 0.22153882 0.11822882

bad_clusters <- which(cluster_prop < 0.01 / k)
if (length(bad_clusters) > 0) {
    message("For k: ", k, " we are dropping small clusters: ", paste(names(bad_clusters), collapse = ", "))
    spe <- spe[, !colData(spe)[[clustering_no_na]] %in% as.integer(names(bad_clusters))]
    updated_t <- find_t(colData(spe)[[clustering_no_na]], 0.01)
    message("initial t: ", initial_t, "; updated t: ", updated_t)
} else{
    updated_t <- initial_t
}

set.seed(7)
fasthplus <- hpb(D = reduced_dim_no_na, L = clustering_no_na, t = updated_t, r = 30)
results <- data.frame(k = k, fasthplus = fasthplus)
write.table(results, file = here::here("processed-data", "08_clustering", "cluster_diagnostics", "fasthplus", "fasthplus_results_nnSVG_precast.csv"),
            append = TRUE,
            row.names = FALSE,
            col.names =!file.exists(here("processed-data", "08_clustering", "cluster_diagnostics", "fasthplus", "fasthplus_results_nnSVG_precast.csv")))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
