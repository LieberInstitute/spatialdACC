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
})

load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))
colnames(spe) <- spe$key

K <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
load(file = here("processed-data", "08_clustering", "PRECAST", paste0("nnSVG_PRECASTObj_",K,".Rdata")))

PRECASTObj <- SelectModel(PRECASTObj)
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

# Merge with spe object
cluster_df <- seuInt@meta.data |>
    mutate(cluster = factor(cluster)) |>
    rename_with(~ paste0("PRECAST_", .x)) |>
    rownames_to_column(var = "key")

col_data_df <- colData(spe) |>
    data.frame() |>
    left_join(cluster_df, by="key")

rownames(col_data_df) <- colnames(spe)
colData(spe)$PRECAST_cluster <- col_data_df$PRECAST_cluster

precast_name <- paste0("nnSVG_PRECAST_captureArea_", K)

cluster_export(
    spe,
    "PRECAST_cluster",
    cluster_dir = here::here("processed-data", "08_clustering", "PRECAST", precast_name)
)

brains <- unique(spe$brnum)

pdf(file = here::here("plots", "08_clustering", "PRECAST", paste0(precast_name, ".pdf")), width = 21, height = 20)

for (i in seq_along(brains)){
    speb <- spe[, which(spe$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    print(length(samples))

    if (length(samples) == 1){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "PRECAST_cluster", colors = cols, point_size = 4, ... = paste0("_", brains[i]))
        grid.arrange(p1, nrow = 1)
    } else if (length(samples) == 2){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "PRECAST_cluster", colors = cols, point_size = 3, ... = paste0("_", brains[i]))
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "PRECAST_cluster", colors = cols, point_size = 3, ... = paste0("_", brains[i]))
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "PRECAST_cluster", colors = cols, point_size = 3, ... = paste0("_", brains[i]))
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "PRECAST_cluster", colors = cols, point_size = 3, ... = paste0("_", brains[i]))
        p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "PRECAST_cluster", colors = cols, point_size = 3, ... = paste0("_", brains[i]))
        grid.arrange(p1, p2, p3, nrow = 2)
    }
}

samples <- unique(colData(spe)[, c("sample_id", "brnum")])
rownames(samples) <- NULL

for (i in 1:nrow(samples)) {
    p <- vis_clus(
        spe = spe,
        sampleid = samples$sample_id[i],
        clustervar = "PRECAST_cluster",
        colors = c("FALSE" = "yellow", "TRUE" = "blue"),
        point_size = 2,
        ... = paste0("_", samples$brnum[i])
    )

    p1 <- plotVisium(spe[, which(spe$sample_id == samples$sample_id[i])], spots = FALSE)

    grid.arrange(p, p1, nrow = 1)
}

dev.off()
