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
})

load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))

read_barcoded_csv <- function(x) {
    df <- read.csv(x)
    colnames(df) <- tolower(colnames(df))

    if (colnames(df)[2] == "cluster") {
        colnames(df)[2] <-
            gsub("gene_expression_", "", basename(dirname(x)))
    }
    return(df)
}

PRECAST_import <- function(spe, cluster_dir = file.path(tempdir(), "exported_clusters")) {
    clustering_files <-
        list.files(
            here::here("processed-data", "08_clustering", "PRECAST", precast_name),
            pattern = "clusters.csv",
            all.files = TRUE,
            full.names = TRUE,
            recursive = TRUE
        )
    clusters_list <- lapply(clustering_files, read_barcoded_csv)
    clusters <- Reduce(function(...) merge(..., by = "key", all = TRUE), clusters_list)
    cluster_cols <- which(colnames(clusters) != "key")
    colnames(clusters)[cluster_cols] <- paste0("", colnames(clusters)[cluster_cols])

    colData(spe)[,precast_name] <- clusters[,2]
    return(spe)
}

#import cluster columns into single spe
num_clusters <- c(5:20)
for (k in num_clusters) {

    #import bayesspace harmony clusters
    bayesSpace_name <- paste0("bayesSpace_captureArea_", k)
    spe <- cluster_import(
        spe,
        cluster_dir = here::here("processed-data", "08_clustering", "BayesSpace", "preprocess_harmony", bayesSpace_name),
        prefix = "harmony_"
    )

    #import bayesspace mnn clusters
    spe <- cluster_import(
        spe,
        cluster_dir = here::here("processed-data", "08_clustering", "BayesSpace", "preprocess_mnn", bayesSpace_name),
        prefix = "mnn_"
    )

    #import precast clusters
    #i re wrote the key when using PRECAST, just add the cluster columns manually
    precast_name <- paste0("PRECAST_captureArea_", k)
    spe <- PRECAST_import(
        spe,
        cluster_dir = here::here("processed-data", "08_clustering", "PRECAST", precast_name)
    )
}

#find 13 spots that PRECAST filtered out
non_na_indices <- !is.na(colData(spe)$PRECAST_captureArea_5)

clustering_columns <- colData(spe)[,c(49:96)]
column_names <- colnames(colData(spe))[c(49:96)]

hplus.list <- list()

for (i in seq_along(clustering_columns)) {
    current_clustering <- clustering_columns[, i]
    current_colname <- column_names[i]

    # Remove NAs from the current_clustering
    current_clustering_no_na <- current_clustering[non_na_indices]
    reduced_dim_no_na <- reducedDim(spe, "pp-GLM-PCA")[non_na_indices,]

    #calculate HPE using bootstrap for current clustering
    hplus <- hpb(reduced_dim_no_na, current_clustering_no_na, t = 0.03*length(current_clustering_no_na))

    #add HPE of current clustering to list, named under current clustering
    hplus.list[[current_colname]] <- hplus
}

save(hplus.list,file=here::here("plots","08_clustering","cluster_diagnostics","discordance.rda"))

