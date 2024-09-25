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

spe$logcounts_RORB <- logcounts(spe)[which(rowData(spe)$gene_name=="RORB"),]

spe.temp <- spe
brains <- unique(spe.temp$brnum)

pdf(file = here::here("plots", "18_PsychENCODE_NMF", "RORB_spot_plots_dACC.pdf"),
    width = 21, height = 20)

for (j in seq_along(brains)){
    speb <- spe.temp[, which(spe.temp$brnum == brains[j])]
    samples <- unique(speb$sample_id)
    print(length(samples))

    if (length(samples) == 1){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="logcounts_RORB", point_size = 9) |> add_ground(var="layer", stroke=0.5, point_size = 9) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        grid.arrange(p1, nrow = 1)
    } else if (length(samples) == 2){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))

        spe_3 <- speb[, which(speb$sample_id == samples[3])]
        p3 <- make_escheR(spe_3) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[3]))
        grid.arrange(p1, p2, p3, nrow = 2)
    }
}

dev.off()


# load data
spe <- spatialLIBD::fetch_data(type = "spe")

spe$logcounts_RORB <- logcounts(spe)[which(rowData(spe)$gene_name=="RORB"),]

spe.temp <- spe
brains <- unique(spe.temp$subject)

spe.temp$layer_guess_reordered <- unfactor(spe.temp$layer_guess_reordered)
spe.temp$layer_guess_reordered[spe.temp$layer_guess_reordered == "Layer1"] <- "L1"
spe.temp$layer_guess_reordered[spe.temp$layer_guess_reordered == "Layer2"] <- "L2"
spe.temp$layer_guess_reordered[spe.temp$layer_guess_reordered == "Layer3"] <- "L3"
spe.temp$layer_guess_reordered[spe.temp$layer_guess_reordered == "Layer4"] <- "L4"
spe.temp$layer_guess_reordered[spe.temp$layer_guess_reordered == "Layer5"] <- "L5"
spe.temp$layer_guess_reordered[spe.temp$layer_guess_reordered == "Layer6"] <- "L6"
# remove NA
spe.temp <- spe.temp[,!is.na(spe.temp$layer_guess_reordered)]

spe.temp$layer <- spe.temp$layer_guess_reordered

pdf(file = here::here("plots", "18_PsychENCODE_NMF", "RORB_spot_plots_DLPFC_12.pdf"),
    width = 21, height = 20)

for (j in seq_along(brains)){
    speb <- spe.temp[, which(spe.temp$subject == brains[j])]
    samples <- unique(speb$sample_id)
    print(length(samples))

    if (length(samples) == 1){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="logcounts_RORB", point_size = 9) |> add_ground(var="layer", stroke=0.5, point_size = 9) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        grid.arrange(p1, nrow = 1)
    } else if (length(samples) == 2){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))

        spe_3 <- speb[, which(speb$sample_id == samples[3])]
        p3 <- make_escheR(spe_3) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[3]))
        grid.arrange(p1, p2, p3, nrow = 2)
    } else if (length(samples) == 4){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))

        spe_3 <- speb[, which(speb$sample_id == samples[3])]
        p3 <- make_escheR(spe_3) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[3]))

        spe_4 <- speb[, which(speb$sample_id == samples[4])]
        p4 <- make_escheR(spe_4) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[3]))
        grid.arrange(p1, p2, p3, nrow = 2)
    }
}

dev.off()


# load data
spe <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium")

spe$logcounts_RORB <- logcounts(spe)[which(rowData(spe)$gene_name=="RORB"),]

spe.temp <- spe
brains <- unique(spe.temp$subject)

# create spatial labels for DLPFC_30
spe.temp$BayesSpace_harmony_09[spe.temp$BayesSpace_harmony_09 == 3] <- "L2"
spe.temp$BayesSpace_harmony_09[spe.temp$BayesSpace_harmony_09 == 8] <- "L4"
spe.temp$BayesSpace_harmony_09[spe.temp$BayesSpace_harmony_09 == 7] <- "L6"
spe.temp$BayesSpace_harmony_09[spe.temp$BayesSpace_harmony_09 == 5] <- "L3"
spe.temp$BayesSpace_harmony_09[spe.temp$BayesSpace_harmony_09 == 6] <- "WM"
spe.temp$BayesSpace_harmony_09[spe.temp$BayesSpace_harmony_09 == 4] <- "L5"
spe.temp$BayesSpace_harmony_09[spe.temp$BayesSpace_harmony_09 == 2] <- "L1"
spe.temp$BayesSpace_harmony_09[spe.temp$BayesSpace_harmony_09 == 1] <- "meninges"
spe.temp$BayesSpace_harmony_09[spe.temp$BayesSpace_harmony_09 == 9] <- "WM"
# remove meninges
spe.temp <- spe.temp[ , which(spe.temp$BayesSpace_harmony_09 != "meninges")]

spe.temp$layer <- spe.temp$BayesSpace_harmony_09

pdf(file = here::here("plots", "18_PsychENCODE_NMF", "RORB_spot_plots_DLPFC_30.pdf"),
    width = 21, height = 20)

for (j in seq_along(brains)){
    speb <- spe.temp[, which(spe.temp$subject == brains[j])]
    samples <- unique(speb$sample_id)
    print(length(samples))

    if (length(samples) == 1){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="logcounts_RORB", point_size = 9) |> add_ground(var="layer", stroke=0.5, point_size = 9) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        grid.arrange(p1, nrow = 1)
    } else if (length(samples) == 2){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))

        spe_3 <- speb[, which(speb$sample_id == samples[3])]
        p3 <- make_escheR(spe_3) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[3]))
        grid.arrange(p1, p2, p3, nrow = 2)
    } else if (length(samples) == 4){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))

        spe_3 <- speb[, which(speb$sample_id == samples[3])]
        p3 <- make_escheR(spe_3) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[3]))

        spe_4 <- speb[, which(speb$sample_id == samples[4])]
        p4 <- make_escheR(spe_4) |> add_fill(var="logcounts_RORB", point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[3]))
        grid.arrange(p1, p2, p3, nrow = 2)
    }
}

dev.off()
