setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(escheR)
library(here)
library(gridExtra)

# load data
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata"))

# spot plots for VAT1L, ADRA1A , GABRQ, ITGA4, POU3F1, BMP3

spe$counts_BMP3 <- counts(spe)[which(rowData(spe)$gene_name=="BMP3"),]

spe.temp <- spe
brains <- unique(spe.temp$brnum)

pdf(file = here::here("plots", "16_VENs_analysis", "BMP3_spot_plots_clusters.pdf"),
    width = 21, height = 20)

for (j in seq_along(brains)){
    speb <- spe.temp[, which(spe.temp$brnum == brains[j])]
    samples <- unique(speb$sample_id)
    print(length(samples))

    if (length(samples) == 1){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="counts_BMP3", point_size = 9) |> add_ground(var="layer", stroke=0.5, point_size = 9) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        grid.arrange(p1, nrow = 1)
    } else if (length(samples) == 2){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="counts_BMP3", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="counts_BMP3", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="counts_BMP3", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="counts_BMP3", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))

        spe_3 <- speb[, which(speb$sample_id == samples[3])]
        p3 <- make_escheR(spe_3) |> add_fill(var="counts_BMP3", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[3]))
        grid.arrange(p1, p2, p3, nrow = 2)
    }
}

dev.off()

