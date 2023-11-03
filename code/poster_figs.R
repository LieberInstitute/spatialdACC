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
               colors = c("#000000", "#009292", "#FF6DB6", "#490092", "#006DDB", "#924900", "#24FF24", "#FF2A00", "#FFFF00"),
               point_size = 2, spatial = FALSE) +
    theme(legend.position="none") +
    ggtitle("PRECAST Clusters") +
    theme(plot.title = element_text(size=15, face="bold"))

ggsave(here("plots", "poster_figs", "fig2a.png"),
       p1

)

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
    point_size = 1.5,
    percentVar = metadata(spe_pseudo)$PCA_var_explained) +
    scale_color_manual(values=c("#000000", "#009292", "#FF6DB6", "#490092", "#006DDB", "#924900", "#24FF24", "#FF2A00", "#FFFF00")) +
    guides(color=guide_legend("cluster"))  +
    ggtitle("PC Scores") +
    theme(plot.title = element_text(size=15)) +
    theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))


vars <- getVarianceExplained(spe_pseudo,
                             variables = c("cluster","sample_id", "sum_sample", "detected_sample")
)


p3 <- plotExplanatoryVariables(vars)  +
    ggtitle("Var. Explained") +
    theme(plot.title = element_text(size=15)) +
    guides(color=guide_legend(""))  +
    scale_color_manual(labels = c("cluster", "sample", "sum", "detected"), values = c("#631879", "#008280", "#808180", "#FDAF91")) +
    theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

wrap_plots((p2 | p3)) &
    theme(plot.tag = element_text(color = "black", size = 20, face="bold"))

ggsave(here("plots", "poster_figs", "fig2b.png"), p2)
ggsave(here("plots", "poster_figs", "fig2c.png"), p3)

#figure 3
colData(spe)[nnSVG_precast_name] <- as.factor(colData(spe)[,nnSVG_precast_name])
colData(spe)[nnSVG_precast_name] <- paste0("clust", colData(spe)[,nnSVG_precast_name])

spe_pseudo <-
    registration_pseudobulk(spe,
                            var_registration = nnSVG_precast_name,
                            var_sample_id = "sample_id",
                            min_ncells = 10
    )

load(
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_DE", paste0(nnSVG_precast_name,".Rdata"))
)
colData(spe)$spatialLIBD <- colData(spe)$registration_variable

sig_genes <- sig_genes_extract(
    n = 20,
    modeling_results = modeling_results,
    model_type = "enrichment",
    sce_layer = spe_pseudo
)

indices <- c(87, 165, 21)

png(here("plots", "poster_figs", "fig3a.png"))

p1 <- layer_boxplot(
    87,
    sig_genes = sig_genes,
    short_title = TRUE,
    sce_layer = spe_pseudo,
    col_bkg_box = "grey80",
    col_bkg_point = "grey40",
    col_low_box = "violet",
    col_low_point = "darkviolet",
    col_high_box = "skyblue",
    col_high_point = "dodgerblue4",
    cex = 2,
    group_var = nnSVG_precast_name,
    assayname = "logcounts"
)

dev.off()

png(here("plots", "poster_figs", "fig3b.png"))

p2 <- layer_boxplot(
    165,
    sig_genes = sig_genes,
    short_title = TRUE,
    sce_layer = spe_pseudo,
    col_bkg_box = "grey80",
    col_bkg_point = "grey40",
    col_low_box = "violet",
    col_low_point = "darkviolet",
    col_high_box = "skyblue",
    col_high_point = "dodgerblue4",
    cex = 2,
    group_var = nnSVG_precast_name,
    assayname = "logcounts"
)

dev.off()

png(here("plots", "poster_figs", "fig3c.png"))


p3 <- layer_boxplot(
    64,
    sig_genes = sig_genes,
    short_title = TRUE,
    sce_layer = spe_pseudo,
    col_bkg_box = "grey80",
    col_bkg_point = "grey40",
    col_low_box = "violet",
    col_low_point = "darkviolet",
    col_high_box = "skyblue",
    col_high_point = "dodgerblue4",
    cex = 2,
    group_var = nnSVG_precast_name,
    assayname = "logcounts"
)

dev.off()


png(here("plots", "poster_figs", "fig3d.png"))

vis_gene(
    spe = spe,
    sampleid = "V12N28-331_D1",
    geneid = rownames(spe)[which(rowData(spe)$gene_name == "RELN")]
)

dev.off()

png(here("plots", "poster_figs", "fig3e.png"))

vis_gene(
    spe = spe,
    sampleid = "V12N28-331_D1",
    geneid = rownames(spe)[which(rowData(spe)$gene_name == "TRABD2A")]
)

dev.off()

png(here("plots", "poster_figs", "fig3f.png"))

vis_gene(
    spe = spe,
    sampleid = "V12N28-331_D1",
    geneid = rownames(spe)[which(rowData(spe)$gene_name == "KRT17")]
)

dev.off()


#fig 4

k=7
nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", k)
spe <- nnSVG_PRECAST_import(
    spe,
    cluster_dir = here::here("processed-data", "08_clustering", "PRECAST", nnSVG_precast_name)
)


colData(spe)[nnSVG_precast_name] <- as.factor(colData(spe)[,nnSVG_precast_name])
colData(spe)[nnSVG_precast_name] <- paste0("clust", colData(spe)[,nnSVG_precast_name])

spe_pseudo <-
    registration_pseudobulk(spe,
                            var_registration = nnSVG_precast_name,
                            var_sample_id = "sample_id",
                            min_ncells = 10
    )

load(
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_DE", paste0(nnSVG_precast_name,".Rdata"))
)
colData(spe)$spatialLIBD <- colData(spe)$registration_variable

sig_genes <- sig_genes_extract(
    n = 20,
    modeling_results = modeling_results,
    model_type = "enrichment",
    sce_layer = spe_pseudo
)

png(here("plots", "poster_figs", "fig4a.png"))

p1 <- layer_boxplot(
    53,
    sig_genes = sig_genes,
    short_title = TRUE,
    sce_layer = spe_pseudo,
    col_bkg_box = "grey80",
    col_bkg_point = "grey40",
    col_low_box = "violet",
    col_low_point = "darkviolet",
    col_high_box = "skyblue",
    col_high_point = "dodgerblue4",
    cex = 2,
    group_var = nnSVG_precast_name,
    assayname = "logcounts"
)

dev.off()

k=10
nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", k)
spe <- nnSVG_PRECAST_import(
    spe,
    cluster_dir = here::here("processed-data", "08_clustering", "PRECAST", nnSVG_precast_name)
)


colData(spe)[nnSVG_precast_name] <- as.factor(colData(spe)[,nnSVG_precast_name])
colData(spe)[nnSVG_precast_name] <- paste0("clust", colData(spe)[,nnSVG_precast_name])

spe_pseudo <-
    registration_pseudobulk(spe,
                            var_registration = nnSVG_precast_name,
                            var_sample_id = "sample_id",
                            min_ncells = 10
    )

load(
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_DE", paste0(nnSVG_precast_name,".Rdata"))
)
colData(spe)$spatialLIBD <- colData(spe)$registration_variable

sig_genes <- sig_genes_extract(
    n = 20,
    modeling_results = modeling_results,
    model_type = "enrichment",
    sce_layer = spe_pseudo
)

png(here("plots", "poster_figs", "fig4b.png"))

p1 <- layer_boxplot(
    96,
    sig_genes = sig_genes,
    short_title = TRUE,
    sce_layer = spe_pseudo,
    col_bkg_box = "grey80",
    col_bkg_point = "grey40",
    col_low_box = "violet",
    col_low_point = "darkviolet",
    col_high_box = "skyblue",
    col_high_point = "dodgerblue4",
    cex = 2,
    group_var = nnSVG_precast_name,
    assayname = "logcounts"
)

dev.off()

p1 <- vis_clus(spe, sampleid = "V12N28-331_D1", clustervar = paste0("nnSVG_PRECAST_captureArea_", 10),
               colors = c("#000000", "#009292", "#FF6DB6", "#490092", "#006DDB", "#924900", "#24FF24", "#FF2A00", "#FFFF00", "#00FFFF"),
               point_size = 2, spatial = FALSE) +
    ggtitle("PRECAST Clusters") +
    theme(plot.title = element_text(size=15, face="bold"))


ggsave(here("plots", "poster_figs", "fig4c.png"),
       p1
)



p1 <- vis_clus(spe, sampleid = "V12N28-331_D1", clustervar = paste0("nnSVG_PRECAST_captureArea_", 7),
               colors = c("#000000", "#009292", "#FF6DB6", "#490092", "#006DDB", "#924900", "#24FF24"),
               point_size = 2, spatial = FALSE) +
    ggtitle("PRECAST Clusters") +
    theme(plot.title = element_text(size=15, face="bold"))

ggsave(here("plots", "poster_figs", "fig4d.png"),
       p1
)
