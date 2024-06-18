setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(escheR)
library(here)
library(gridExtra)

# load data
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9.Rdata"))

spe$counts_VAT1L <- counts(spe)[which(rowData(spe)$gene_name=="VAT1L"),]

spe$cluster <- unfactor(spe$cluster)
spe$cluster[spe$cluster == 3] <- "WM1"
spe$cluster[spe$cluster == 8] <- "WM2"
spe$cluster[spe$cluster == 7] <- "WM-CC"
spe$cluster[spe$cluster == 5] <- "L6b"
spe$cluster[spe$cluster == 6] <- "L6a"
spe$cluster[spe$cluster == 4] <- "L5"
spe$cluster[spe$cluster == 2] <- "L3"
spe$cluster[spe$cluster == 1] <- "L2"
spe$cluster[spe$cluster == 9] <- "L1"

spe.temp <- spe
brains <- unique(spe.temp$subject)

pdf(file = here::here("plots", "16_VENs_analysis", "VAT1L_spot_plots_clusters.pdf"),
    width = 21, height = 20)

for (j in seq_along(brains)){
    speb <- spe.temp[, which(spe.temp$subject == brains[j])]
    samples <- unique(speb$sample_id)
    print(length(samples))

    if (length(samples) == 1){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="counts_VAT1L", point_size = 9) |> add_ground(var="cluster", stroke=0.5, point_size = 9) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        grid.arrange(p1, nrow = 1)
    } else if (length(samples) == 2){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="counts_VAT1L", point_size = 4) |> add_ground(var="cluster", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="counts_VAT1L", point_size = 4) |> add_ground(var="cluster", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="counts_VAT1L", point_size = 4) |> add_ground(var="cluster", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="counts_VAT1L", point_size = 4) |> add_ground(var="cluster", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))

        spe_3 <- speb[, which(speb$sample_id == samples[3])]
        p3 <- make_escheR(spe_3) |> add_fill(var="counts_VAT1L", point_size = 4) |> add_ground(var="cluster", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[3]))
        grid.arrange(p1, p2, p3, nrow = 2)
    }
}

dev.off()


# create violin plots to show the distribution of VAT1L expression in each spatial domain
spe_DLPFC_30 <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium")

spe_DLPFC_30$counts_VAT1L <- counts(spe_DLPFC_30)[which(rowData(spe_DLPFC_30)$gene_name=="VAT1L"),]

# create spatial labels for DLPFC_30
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 3] <- "L2"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 8] <- "L4"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 7] <- "L6"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 5] <- "L3"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 6] <- "WM"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 4] <- "L5"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 2] <- "L1"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 1] <- "meninges"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 9] <- "WM"

# remove meninges
spe_DLPFC_30 <- spe_DLPFC_30[,which(spe_DLPFC_30$BayesSpace_harmony_09 != "meninges")]

spe_DLPFC_30$cluster <- spe_DLPFC_30$BayesSpace_harmony_09

spe.temp <- spe_DLPFC_30
brains <- unique(spe.temp$subject)

pdf(file = here::here("plots", "16_VENs_analysis", "VAT1L_spot_plots_clusters_DLPFC_30.pdf"),
    width = 21, height = 20)

for (j in seq_along(brains)){
    speb <- spe.temp[, which(spe.temp$subject == brains[j])]
    samples <- unique(speb$sample_id)
    print(length(samples))

    if (length(samples) == 1){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="counts_VAT1L", point_size = 9) |> add_ground(var="cluster", stroke=0.5, point_size = 9) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        grid.arrange(p1, nrow = 1)
    } else if (length(samples) == 2){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="counts_VAT1L", point_size = 4) |> add_ground(var="cluster", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="counts_VAT1L", point_size = 4) |> add_ground(var="cluster", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="counts_VAT1L", point_size = 4) |> add_ground(var="cluster", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="counts_VAT1L", point_size = 4) |> add_ground(var="cluster", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))

        spe_3 <- speb[, which(speb$sample_id == samples[3])]
        p3 <- make_escheR(spe_3) |> add_fill(var="counts_VAT1L", point_size = 4) |> add_ground(var="cluster", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[3]))
        grid.arrange(p1, p2, p3, nrow = 2)
    }
}

dev.off()
