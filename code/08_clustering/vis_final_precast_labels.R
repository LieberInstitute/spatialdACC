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

load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata"))

brains <- unique(spe$brnum)

pdf(file = here::here("plots", "08_clustering", "PRECAST", "final_clusters.pdf"), width = 21, height = 20)

for (i in seq_along(brains)){
    speb <- spe[, which(spe$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    print(length(samples))

    if (length(samples) == 1){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "layer", colors = cols, spatial = FALSE, point_size = 8, ... = paste0("_", brains[i])) +
            scale_fill_manual(values = c(
                "L2" = "#CD661D",
                "L3" = "#00FA9A",
                "L5" = "#00CED1",
                "L6a" = "#000080",
                "L6b" = "#EED5B7",
                "WM" = "black",
                "L1" = "#CDB7B5"
            ))
        grid.arrange(p1, nrow = 1)
    } else if (length(samples) == 2){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "layer", colors = cols, spatial = FALSE, point_size = 4, ... = paste0("_", brains[i])) +
            scale_fill_manual(values = c(
                "L2" = "#CD661D",
                "L3" = "#00FA9A",
                "L5" = "#00CED1",
                "L6a" = "#000080",
                "L6b" = "#EED5B7",
                "WM" = "black",
                "L1" = "#CDB7B5"
            ))
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "layer", colors = cols, spatial = FALSE, point_size = 4, ... = paste0("_", brains[i])) +
            scale_fill_manual(values = c(
                "L2" = "#CD661D",
                "L3" = "#00FA9A",
                "L5" = "#00CED1",
                "L6a" = "#000080",
                "L6b" = "#EED5B7",
                "WM" = "black",
                "L1" = "#CDB7B5"
            ))
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "layer", colors = cols, spatial = FALSE, point_size = 3, ... = paste0("_", brains[i])) +
            scale_fill_manual(values = c(
                "L2" = "#CD661D",
                "L3" = "#00FA9A",
                "L5" = "#00CED1",
                "L6a" = "#000080",
                "L6b" = "#EED5B7",
                "WM" = "black",
                "L1" = "#CDB7B5"
            ))
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "layer", colors = cols, spatial = FALSE, point_size = 3, ... = paste0("_", brains[i])) +
            scale_fill_manual(values = c(
                "L2" = "#CD661D",
                "L3" = "#00FA9A",
                "L5" = "#00CED1",
                "L6a" = "#000080",
                "L6b" = "#EED5B7",
                "WM" = "black",
                "L1" = "#CDB7B5"
            ))
        p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "layer", colors = cols, spatial = FALSE, point_size = 3, ... = paste0("_", brains[i])) +
            scale_fill_manual(values = c(
                "L2" = "#CD661D",
                "L3" = "#00FA9A",
                "L5" = "#00CED1",
                "L6a" = "#000080",
                "L6b" = "#EED5B7",
                "WM" = "black",
                "L1" = "#CDB7B5"
            ))
        grid.arrange(p1, p2, p3, nrow = 2)
    }
}

samples <- unique(colData(spe)[, c("sample_id", "brnum")])
rownames(samples) <- NULL

for (i in 1:nrow(samples)) {
    p <- vis_clus(
        spe = spe,
        sampleid = samples$sample_id[i],
        clustervar = "layer",
        colors = c("FALSE" = "yellow", "TRUE" = "blue"),
        spatial = FALSE,
        point_size = 4,
        ... = paste0("_", samples$brnum[i])
    )

    p1 <- plotVisium(spe[, which(spe$sample_id == samples$sample_id[i])], spots = FALSE)

    grid.arrange(p, p1, nrow = 1)
}

dev.off()
