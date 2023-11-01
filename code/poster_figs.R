setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("scater")
    library("spatialLIBD")
    library("dplyr")
    library("patchwork")
    library("ggplot2")
    library("purrr")
    library("tidyverse")
    library("gridExtra")
    library("ggspavis")
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


#figure 1 made by madeline valentine

#figure 2

load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))

k=9
nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", k)
spe <- nnSVG_PRECAST_import(
    spe,
    cluster_dir = here::here("processed-data", "08_clustering", "PRECAST", nnSVG_precast_name)
)

p1 <- vis_clus(spe, sampleid = "V12N28-331_D1", clustervar = nnSVG_precast_name,
               colors = c("#c2cea5", "#e9c891", "#c891e9", "#91e9c8", "#f57b9d", "#c6c506", "#c506c6", "#06c6c5", "#f8b200"),
               point_size = 4, spatial = FALSE) +
    theme(legend.position="none") +
    ggtitle("PRECAST Clusters") +
    theme(plot.title = element_text(size=15, face="bold"))

load(
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_pseudobulk", paste0(nnSVG_precast_name,".Rdata"))
)

sum_by_sample <- setNames(aggregate(sum ~ sample_id, colData(spe), sum), c("sample_id", "sum_sample"))
detected_by_sample <- setNames(aggregate(detected ~ sample_id, colData(spe), sum), c("sample_id", "detected_sample"))

col_data_df <- as.data.frame(colData(spe_pseudo))
col_data_df <- left_join(col_data_df, detected_by_sample, by = "sample_id")
col_data_df <- left_join(col_data_df, sum_by_sample, by = "sample_id")

colData(spe_pseudo) <- DataFrame(col_data_df)

colData(spe_pseudo)["cluster"] <- colData(spe_pseudo)[,nnSVG_precast_name]

p2 <- plotPCA(
    spe_pseudo,
    colour_by = "cluster",
    ncomponents = 2,
    point_size = 3,
    percentVar = metadata(spe_pseudo)$PCA_var_explained) +
    scale_color_manual(values=c("#c2cea5", "#e9c891", "#c891e9", "#91e9c8", "#f57b9d", "#c6c506", "#c506c6", "#06c6c5", "#f8b200")) +
    guides(color=guide_legend("cluster"))  +
    ggtitle("PC Scores") +
    theme(plot.title = element_text(size=15))


vars <- getVarianceExplained(spe_pseudo,
                             variables = c("cluster","sample_id", "sum_sample", "detected_sample")
)


p3 <- plotExplanatoryVariables(vars)  +
    ggtitle("PC Variance Explained") +
    theme(plot.title = element_text(size=15))


ggsave(here("plots", "11_differential_expression", "pseudobulk", "nnSVG_precast_pseudobulk", paste0("pseudobulk_PC_",nnSVG_precast_name,".pdf")),
                wrap_plots(p1,p2,p3, nrow=2) + plot_annotation(tag_levels = 'A'))
