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

load(file = here("processed-data", "18_PsychENCODE_NMF", "spe_DLPFC_30_NMF.Rdata"))
load(file = here("processed-data", "18_PsychENCODE_NMF", "spe_DLPFC_12_NMF.Rdata"))
load(file = here("processed-data", "18_PsychENCODE_NMF", "spe_dACC_NMF.Rdata"))

spe_DLPFC_12$layer_guess_reordered <- unfactor(spe_DLPFC_12$layer_guess_reordered)
spe_DLPFC_12$layer_guess_reordered[spe_DLPFC_12$layer_guess_reordered == "Layer1"] <- "L1"
spe_DLPFC_12$layer_guess_reordered[spe_DLPFC_12$layer_guess_reordered == "Layer2"] <- "L2"
spe_DLPFC_12$layer_guess_reordered[spe_DLPFC_12$layer_guess_reordered == "Layer3"] <- "L3"
spe_DLPFC_12$layer_guess_reordered[spe_DLPFC_12$layer_guess_reordered == "Layer4"] <- "L4"
spe_DLPFC_12$layer_guess_reordered[spe_DLPFC_12$layer_guess_reordered == "Layer5"] <- "L5"
spe_DLPFC_12$layer_guess_reordered[spe_DLPFC_12$layer_guess_reordered == "Layer6"] <- "L6"

#remove NA values
spe_DLPFC_12 <- spe_DLPFC_12[,!is.na(spe_DLPFC_12$layer_guess_reordered)]

# rename to "layer"
colData(spe_DLPFC_12)$layer <- colData(spe_DLPFC_12)$layer_guess_reordered

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

# remove meninges from DLPFC_30
spe_DLPFC_30 <- spe_DLPFC_30[,spe_DLPFC_30$BayesSpace_harmony_09 != "meninges"]

DLPFC_30_layer <- "BayesSpace_harmony_09"
DLPFC_12_layer <- "layer"
dACC_layer <- "layer"

# create a function to make a list of 3 plots
# the input is NMF factor, sample id for DLPFC_12, sample id for DLPFC 30, and sample id for dACC
# the first plot is escher plot for DLPFC_12 with add_ground(var = DLPFC_12_layer)
# the second plot is escher plot for DLPFC_30 with add_ground(var = DLPFC_30_layer)
# the third plot is escher plot for dACC with add_ground(var = dACC_layer)
# the fill for the escher plots is the NMF factor

make_plots <- function(factor, DLPFC_12_id, DLPFC_30_id, dACC_id){

    spe_1 <- spe_DLPFC_12[, which(spe_DLPFC_12$sample_id == DLPFC_12_id)]
    p1 <- make_escheR(spe_1) |>
        add_fill(var=factor, point_size = 1) |>
        add_ground(var=DLPFC_12_layer, stroke=0.15, point_size = 1) +
        scale_color_manual(values = c(
            "L2" = "#E41A1C",   # Bright red
            "L3" = "#377EB8",   # Strong blue
            "L5" = "#4DAF4A",   # Vivid green
            "L4" = "#984EA3",  # Purple
            "L6" = "#FF7F00",  # Orange
            "WM" = "#FFFF99",  # Light yellow
            "L1" = "#00CED1"    # Dark turquoise
        )) +
        scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", DLPFC_12_id)) +
        theme(
            legend.text = element_text(size = 5), # Smaller text size
            legend.title = element_text(size = 8), # Smaller title size
            legend.key.size = unit(0.5, "lines"), # Smaller key size
            legend.spacing.x = unit(0.3, 'cm'), # Adjust spacing between keys
            legend.spacing.y = unit(0.3, 'cm'))  # Adjust spacing between rows

    spe_2 <- spe_DLPFC_30[, which(spe_DLPFC_30$sample_id == DLPFC_30_id)]
    p2 <- make_escheR(spe_2) |>
        add_fill(var=factor, point_size = 1) |>
        add_ground(var=DLPFC_30_layer, stroke=0.15, point_size = 1) +
        scale_color_manual(values = c(
            "L2" = "#E41A1C",   # Bright red
            "L3" = "#377EB8",   # Strong blue
            "L5" = "#4DAF4A",   # Vivid green
            "L4" = "#984EA3",  # Purple
            "L6" = "#FF7F00",  # Orange
            "WM" = "#FFFF99",  # Light yellow
            "L1" = "#00CED1"    # Dark turquoise
        )) +
        scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", DLPFC_30_id)) +
        labs(color = "layer") +
        theme(
            legend.text = element_text(size = 5), # Smaller text size
            legend.title = element_text(size = 8), # Smaller title size
            legend.key.size = unit(0.5, "lines"), # Smaller key size
            legend.spacing.x = unit(0.3, 'cm'), # Adjust spacing between keys
            legend.spacing.y = unit(0.3, 'cm'))  # Adjust spacing between rows


    spe_3 <- spe_dACC[, which(spe_dACC$sample_id == dACC_id)]
    p3 <- make_escheR(spe_3) |>
        add_fill(var=factor, point_size = 1) |>
        add_ground(var=dACC_layer, stroke=0.15, point_size = 1) +
        scale_color_manual(values = c(
            "L2" = "#E41A1C",   # Bright red
            "L3" = "#377EB8",   # Strong blue
            "L5" = "#4DAF4A",   # Vivid green
            "L6a" = "#E37100",  # Orange
            "L6b" = "#FFAE00",  # Orange
            "WM" = "#FFFF99",  # Light yellow
            "L1" = "#00CED1"    # Dark turquoise
        )) +
        scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", dACC_id)) +
        labs(color = "layer") +
        theme(
            legend.text = element_text(size = 5), # Smaller text size
            legend.title = element_text(size = 8), # Smaller title size
            legend.key.size = unit(0.5, "lines"), # Smaller key size
            legend.spacing.x = unit(0.3, 'cm'), # Adjust spacing between keys
            legend.spacing.y = unit(0.3, 'cm'))  # Adjust spacing between rows


    return(list(p1, p2, p3))
}


make_plots_9 <- function(factor, DLPFC_12_id, DLPFC_30_id){

    spe_1 <- spe_DLPFC_12[, which(spe_DLPFC_12$sample_id == DLPFC_12_id)]
    p1 <- make_escheR(spe_1) |>
        add_fill(var=factor, point_size = 1) |>
        add_ground(var=DLPFC_12_layer, stroke=0.15, point_size = 2) +
        scale_color_manual(values = c(
            "L2" = "#E41A1C",   # Bright red
            "L3" = "#377EB8",   # Strong blue
            "L5" = "#4DAF4A",   # Vivid green
            "L4" = "#984EA3",  # Purple
            "L6" = "#FF7F00",  # Orange
            "WM" = "#FFFF99",  # Light yellow
            "L1" = "#00CED1"    # Dark turquoise
        )) +
        scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", DLPFC_12_id)) +
        theme(
            legend.text = element_text(size = 5), # Smaller text size
            legend.title = element_text(size = 8), # Smaller title size
            legend.key.size = unit(0.5, "lines"), # Smaller key size
            legend.spacing.x = unit(0.3, 'cm'), # Adjust spacing between keys
            legend.spacing.y = unit(0.3, 'cm'))  # Adjust spacing between rows

    spe_2 <- spe_DLPFC_30[, which(spe_DLPFC_30$sample_id == DLPFC_30_id)]
    p2 <- make_escheR(spe_2) |>
        add_fill(var=factor, point_size = 1) |>
        add_ground(var=DLPFC_30_layer, stroke=0.15, point_size = 2) +
        scale_color_manual(values = c(
            "L2" = "#E41A1C",   # Bright red
            "L3" = "#377EB8",   # Strong blue
            "L5" = "#4DAF4A",   # Vivid green
            "L4" = "#984EA3",  # Purple
            "L6" = "#FF7F00",  # Orange
            "WM" = "#FFFF99",  # Light yellow
            "L1" = "#00CED1"    # Dark turquoise
        )) +
        scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", DLPFC_30_id)) +
        labs(color = "layer") +
        theme(
            legend.text = element_text(size = 5), # Smaller text size
            legend.title = element_text(size = 8), # Smaller title size
            legend.key.size = unit(0.5, "lines"), # Smaller key size
            legend.spacing.x = unit(0.3, 'cm'), # Adjust spacing between keys
            legend.spacing.y = unit(0.3, 'cm'))  # Adjust spacing between rows


    return(list(p1, p2))
}

plots_9 <- make_plots_9("nmf9", "151673", "Br2743_ant")
plots_10 <- make_plots("nmf10", "151507", "Br6522_ant", "V12J03-002_C1")
plots_15 <- make_plots("nmf15", "151509", "Br6432_mid", "V12N28-331_C1")
plots_16 <- make_plots("nmf16", "151671", "Br6471_ant", "V12N28-334_B1")
plots_20 <- make_plots("nmf20", "151674", "Br3942_ant", "V12N28-332_D1")
plots_23 <- make_plots("nmf23", "151508", "Br6471_mid", "V12J03-002_A1")
plots_32 <- make_plots("nmf32", "151675", "Br3942_post", "V12N28-331_D1")
plots_35 <- make_plots("nmf35", "151673", "Br2720_post", "V12Y31-080_C1")
plots_43 <- make_plots("nmf43", "151676", "Br2720_mid", "V12N28-334_D1")

pdf(here::here("plots", "18_PsychENCODE_NMF", "compare_NMF.pdf"), width = 10, height = 10)

grid.arrange(grobs = plots_9, nrow = 1)
grid.arrange(grobs = plots_10, nrow = 1)
grid.arrange(grobs = plots_15, nrow = 1)
grid.arrange(grobs = plots_16, nrow = 1)
grid.arrange(grobs = plots_20, nrow = 1)
grid.arrange(grobs = plots_23, nrow = 1)
grid.arrange(grobs = plots_32, nrow = 1)
grid.arrange(grobs = plots_35, nrow = 1)
grid.arrange(grobs = plots_43, nrow = 1)

dev.off()

