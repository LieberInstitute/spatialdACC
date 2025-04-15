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
    library(bluster)
    library(patchwork)
})

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

Sys.time()

load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))

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

    #import precast nnSVG clusters
    #i re wrote the key when using PRECAST, just add the cluster columns manually
    nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", k)
    spe <- nnSVG_PRECAST_import(
        spe,
        cluster_dir = here::here("processed-data", "08_clustering", "PRECAST", nnSVG_precast_name)
    )
}

Sys.time()

#find 13 spots that PRECAST filtered out
non_na_indices <- !is.na(colData(spe)$PRECAST_captureArea_5)

clustering_columns <- colData(spe)[,c(49:112)]
column_names <- colnames(colData(spe))[c(49:112)]

purity.data.list <- list()

for (i in seq_along(clustering_columns)) {
    current_clustering <- clustering_columns[, i]
    current_colname <- column_names[i]

    # Remove NAs from the current_clustering
    current_clustering_no_na <- current_clustering[non_na_indices]
    reduced_dim_no_na <- reducedDim(spe, "pp-GLM-PCA")[non_na_indices,]

    #create dataframe of purity data for current clustering
    #purity and max calculated for each gene in current clustering
    purity.approx <- neighborPurity(reduced_dim_no_na, clusters = current_clustering_no_na)
    purity.data <- as.data.frame(purity.approx)
    purity.data$maximum <- factor(purity.data$maximum)
    purity.data$cluster <- current_clustering_no_na

    #add dataframe of current clustering to list, named under current clustering
    purity.data.list[[current_colname]] <- purity.data
}

avg_purity <- data.frame(clustering = character(), cluster = integer(), avg_purity = numeric())
for (colname in names(purity.data.list)) {
    purity.data <- purity.data.list[[colname]]
    #for current clustering, compute average purity of genes per cluster
    avg_purity_col <- aggregate(purity.data$purity, by = list(purity.data$cluster), FUN = mean)
    colnames(avg_purity_col) <- c("cluster", "avg_purity")
    avg_purity_col$clustering <- colname
    avg_purity <- rbind(avg_purity, avg_purity_col)
}

avg_purity$clustering <- factor(avg_purity$clustering, levels = names(purity.data.list))
# Add a column for the algorithm
avg_purity$algorithm <- factor(ifelse(grepl("harmony", avg_purity$clustering), "Harmony - BS",
                                      ifelse(grepl("mnn", avg_purity$clustering), "MNN - BS",
                                             ifelse(grepl("nnSVG", avg_purity$clustering), "nnSVG PRECAST",
                                             "PRECAST"))))
#dataframe has x values for clustering with x clusters
#for example, box plot for harmony_bayesSpace_captureArea_5 is created with 5 average purity values for each of the 5 clusters
save(avg_purity,file=here::here("plots","08_clustering","cluster_diagnostics","purity_boxplot.rda"))


boxplot <- ggplot(avg_purity, aes(x = clustering, y = avg_purity,fill=algorithm)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    #facet_wrap(~ algorithm, strip.position = "bottom")+
    scale_fill_manual(values = c("Harmony - BS" = "#f88379", "nnSVG PRECAST" = "#afeeee",
                                 "PRECAST" = "#9A32CD", "MNN - BS" = "#f9e09c"))+
    ggtitle("Average cluster purity for different clusterings")

pdf(here::here('plots', '08_clustering', 'cluster_diagnostics','avg_purity_boxplot_with_nnSVG.pdf'))
print(boxplot)
dev.off()



# just nnSVG PRECAST diagnostics

avg_purity_subset <- avg_purity[which(avg_purity$algorithm == "nnSVG PRECAST"),]

# just nnSVG PRECAST diagnostics

avg_purity_subset <- avg_purity[which(avg_purity$algorithm == "nnSVG PRECAST"),]
avg_purity_subset$clustering <- sub(".*nnSVG_PRECAST_captureArea_", "", avg_purity_subset$clustering)
avg_purity_subset$clustering <- as.numeric(avg_purity_subset$clustering)
avg_purity_subset$clustering <- as.factor(avg_purity_subset$clustering)

boxplot <- ggplot(avg_purity_subset, aes(x = clustering, y = avg_purity)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggtitle("Average cluster purity nnSVG-Guided PRECAST") +
    theme_bw() +
    xlab("Number of Clusters") +
    ylab("Average Cluster Purity")


df_nnSVG_precast <- read.table(file=here("processed-data", "08_clustering", "cluster_diagnostics", "fasthplus", "fasthplus_results_nnSVG_precast.csv"), header=T)
p1 <- ggplot(data = df_nnSVG_precast, aes(x = k, y = fasthplus, group = 1)) +
    geom_line() +
    geom_point() +
    ylab(expression(H^{
        "+"
    })) +
    theme_bw() +
    ggtitle("H+ Discordance nnSVG-Guided PRECAST")

png(here("plots", "08_clustering", "diagnostics_nnSVG_PRECAST.png"), height=8, width = 12, units = "in", res = 300)
wrap_plots(boxplot,p1, nrow=2) + plot_annotation(tag_levels = 'A')
dev.off()

