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

load(file = here("processed-data", "15_cross_region_snRNA-seq", "spe_DLPFC_30_NMF.Rdata"))
load(file = here("processed-data", "15_cross_region_snRNA-seq", "spe_DLPFC_12_NMF.Rdata"))
load(file = here("processed-data", "15_cross_region_snRNA-seq", "spe_dACC_NMF.Rdata"))

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

# load layer data
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9.Rdata"))

spe$PRECAST_cluster <- unfactor(spe$PRECAST_cluster)
spe$PRECAST_cluster[spe$PRECAST_cluster == 3] <- "WM"
spe$PRECAST_cluster[spe$PRECAST_cluster == 8] <- "WM"
spe$PRECAST_cluster[spe$PRECAST_cluster == 7] <- "WM"
spe$PRECAST_cluster[spe$PRECAST_cluster == 5] <- "L6"
spe$PRECAST_cluster[spe$PRECAST_cluster == 6] <- "L6"
spe$PRECAST_cluster[spe$PRECAST_cluster == 4] <- "L5"
spe$PRECAST_cluster[spe$PRECAST_cluster == 2] <- "L3"
spe$PRECAST_cluster[spe$PRECAST_cluster == 1] <- "L2"
spe$PRECAST_cluster[spe$PRECAST_cluster == 9] <- "L1"

spe_dACC$PRECAST_cluster <- spe$PRECAST_cluster

# add each proj column to colData(spe_DLPFC_30)
for (i in 1:75){
    colData(spe_DLPFC_30)[[paste0("NMF_",i)]] <- reducedDims(spe_DLPFC_30)$NMF_proj[,i]
}

# add each proj column to colData(spe_DLPFC_12)
for (i in 1:75){
    colData(spe_DLPFC_12)[[paste0("NMF_",i)]] <- reducedDims(spe_DLPFC_12)$NMF_proj[,i]
}

# add each proj column to colData(spe_dACC)
for (i in 1:75){
    colData(spe_dACC)[[paste0("NMF_",i)]] <- reducedDims(spe_dACC)$NMF_proj[,i]
}

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
dACC_layer <- "PRECAST_cluster"

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
        add_ground(var=DLPFC_12_layer, stroke=0.3, point_size = 1) +
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
        add_ground(var=DLPFC_30_layer, stroke=0.3, point_size = 1) +
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
        add_ground(var=dACC_layer, stroke=0.3, point_size = 1) +
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

# from unscaled heatmap
plots_53 <- make_plots("NMF_53", "151509", "Br2720_ant", "V12N28-331_D1")

# from rowscale
plots_21 <- make_plots("NMF_21", "151509", "Br2720_ant", "V12N28-331_D1")
plots_25 <- make_plots("NMF_25", "151510", "Br2720_ant", "V12N28-331_D1")

# from columnscale
plots_10 <- make_plots("NMF_10", "151673", "Br3942_post", "V12N28-334_C1")
plots_7 <- make_plots("NMF_7", "151674", "Br2743_post", "V12N28-332_D1")
plots_5 <- make_plots("NMF_5", "151675", "Br2743_mid", "V12N28-331_B1")

# try to find other layer patterns that show up in dACC
#65, 44,59
#16, 30
#52, 42, 46, 18, 72, 49
#43, 75
#13, 6, 24
#54, 51, 55

plots_33 <- make_plots("NMF_33", "151669", "Br3942_post", "V12J03-002_A1")

pdf(here::here("plots", "15_cross_region_snRNA-seq", "compare_NMF.pdf"), width = 10, height = 10)

grid.arrange(grobs = plots_5, nrow = 1)
grid.arrange(grobs = plots_33, nrow = 1)
grid.arrange(grobs = plots_53, nrow = 1)
grid.arrange(grobs = plots_21, nrow = 1)
grid.arrange(grobs = plots_25, nrow = 1)
grid.arrange(grobs = plots_10, nrow = 1)
grid.arrange(grobs = plots_7, nrow = 1)

dev.off()

