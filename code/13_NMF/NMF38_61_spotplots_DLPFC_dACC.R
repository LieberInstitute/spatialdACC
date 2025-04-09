setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(scater)
library(RcppML)
library(ggspavis)
library(here)
library(scRNAseq)
library(Matrix)
library(scran)
library(scuttle)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(igraph)
library(bluster)
library(patchwork)
library(cowplot)
library(projectR)
library(spatialLIBD)
library(gridExtra)
library(escheR)
library(dplyr)
library(tidyverse)

load(file = here("processed-data", "13_NMF", "spe_NMF.Rdata"))
spe_dACC <- spe

load(file = here("processed-data", "13_NMF", "spe_NMF_DLPFC_30.Rdata"))
spe_DLPFC_30 <- spe

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
spe_DLPFC_30 <- spe_DLPFC_30[ , which(spe_DLPFC_30$BayesSpace_harmony_09 != "meninges")]

spe_1 <- spe_dACC[, which(spe_dACC$sample_id == "V12Y31-080_B1")]
p1 <- make_escheR(spe_1) |> add_fill(var="nmf38", point_size = 2) |> add_ground(var="layer", stroke=0.3, point_size = 2) +
    scale_color_manual(values = c(
        "L2" = "#377EB8",
        "L3" = "#4DAF4A",
        "L5" = "#FFD700",
        "L6b" = "#c46200",
        "L6a" = "#FFC18A",
        "WM" = "#1A1A1A",
        "L1" = "#F0027F"
    ),guide="none") +
    scale_fill_gradient(low = "white", high = "black", guide="none") +
    ggtitle("dACC")

p2 <- make_escheR(spe_1) |> add_fill(var="nmf61", point_size = 2) |> add_ground(var="layer", stroke=0.3, point_size = 2) +
    scale_color_manual(values = c(
        "L2" = "#377EB8",
        "L3" = "#4DAF4A",
        "L5" = "#FFD700",
        "L6b" = "#c46200",
        "L6a" = "#FFC18A",
        "WM" = "#1A1A1A",
        "L1" = "#F0027F"
    ), guide="none") +
    scale_fill_gradient(low = "white", high = "black", guide="none")

spe_2 <- spe_DLPFC_30[, which(spe_DLPFC_30$sample_id == "Br8667_mid")]
p3 <- make_escheR(spe_2) |> add_fill(var="nmf38", point_size = 2) |> add_ground(var="BayesSpace_harmony_09", stroke=0.3, point_size = 2) +
    scale_color_manual(values = c(
        "L2" = "#377EB8",
        "L3" = "#4DAF4A",
        "L5" = "#FFD700",
        "L4" = "#984EA3",
        "L6" = "#FF7F00",
        "WM" = "#1A1A1A",
        "L1" = "#F0027F"
    ), guide="none") +
    scale_fill_gradientn(colors=c("white","black"),
                         breaks=c(0,0.001101739),labels=c("low","high"),
                         limits=c(0,0.001101739), name="NMF38") +
    ggtitle("dlPFC")

p4 <- make_escheR(spe_2) |> add_fill(var="nmf61", point_size = 2) |> add_ground(var="BayesSpace_harmony_09", stroke=0.3, point_size = 2) +
    scale_color_manual(values = c(
        "L2" = "#377EB8",
        "L3" = "#4DAF4A",
        "L5" = "#FFD700",
        "L4" = "#984EA3",
        "L6" = "#FF7F00",
        "WM" = "#1A1A1A",
        "L1" = "#F0027F"
    ), guide="none") +
    scale_fill_gradientn(colors=c("white","black"),
                         breaks=c(0,0.0233128),labels=c("low","high"),
                         limits=c(0,0.0233128), name="NMF61")

pdf(file = here::here("plots", "13_NMF", "NMF_spotplots_DLPFC_dACC.pdf"), height = 10, width = 10)
wrap_plots(list(p1,p3,p2,p4), nrow = 2)
dev.off()



