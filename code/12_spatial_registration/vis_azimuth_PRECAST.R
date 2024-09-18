setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("PRECAST")
    library("tictoc")
    library("dplyr")
    library("purrr")
    library("tidyverse")
    library("spatialLIBD")
    library("gridExtra")
    library("ggspavis")
    library("SingleCellExperiment")
    library("stringr")
    library("ComplexHeatmap")
    library("cowplot")
})

# load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))
# colnames(spe) <- spe$key
#
# K <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# load(file = here("processed-data", "08_clustering", "PRECAST", paste0("nnSVG_PRECASTObj_",K,".Rdata")))
#
# PRECASTObj <- SelectModel(PRECASTObj)
# seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")
#
# # Merge with spe object
# cluster_df <- seuInt@meta.data |>
#     mutate(cluster = factor(cluster)) |>
#     rename_with(~ paste0("PRECAST_", .x)) |>
#     rownames_to_column(var = "key")
#
# col_data_df <- colData(spe) |>
#     data.frame() |>
#     left_join(cluster_df, by="key")
#
# rownames(col_data_df) <- colnames(spe)
# colData(spe)$PRECAST_cluster <- col_data_df$PRECAST_cluster

# load spe for k=9 without WM-CC
load(
    file = here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata")
)

precast_name <- paste0("nnSVG_PRECAST_captureArea_", 9)

load(here("processed-data", "12_spatial_registration",paste0("azimuth",precast_name,".rds")))

samples <- unique(colData(spe)[, c("sample_id", "brnum")])
rownames(samples) <- NULL

cor_layer <- cor_layer[order(rownames(cor_layer)), order(colnames(cor_layer))]
p1 <- grid.grabExpr(draw(Heatmap(cor_layer, heatmap_legend_param = list(title = ""))))

pdf(file = here("plots","12_spatial_registration","azimuth",paste0("azimuth_",precast_name,"_combined",".pdf")), height=12, width=12)

for (i in 1:nrow(samples)) {

	p <- vis_clus(
        spe = spe,
        sampleid = samples$sample_id[i],
        clustervar = "layer",
	spatial = FALSE,
        point_size = 1,
        ... = paste0("_", samples$brnum[i])
    )  +
    #reduce size of ggtitle
    theme(plot.title = element_text(size = 12)) +
	#increase legend point size
    guides(color = guide_legend(override.aes = list(size = 3)))

    print(plot_grid(p, p1, labels = "AUTO", ncol = 2))
}

dev.off()
