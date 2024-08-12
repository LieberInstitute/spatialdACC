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

# load dACC data
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9.Rdata"))
spe_clusters <- spe
load(file = here("processed-data", "13_NMF", "spe_NMF.Rdata"))


# add each proj column to colData(spe)
for (i in 1:100){
    colData(spe)[[paste0("NMF_",i)]] <- reducedDims(spe)$NMF_proj[,i]
}

spe_clusters$PRECAST_cluster <- unfactor(spe_clusters$PRECAST_cluster)
spe_clusters$PRECAST_cluster[spe_clusters$PRECAST_cluster == 3] <- "WM1"
spe_clusters$PRECAST_cluster[spe_clusters$PRECAST_cluster == 8] <- "WM2"
spe_clusters$PRECAST_cluster[spe_clusters$PRECAST_cluster == 7] <- "WM-CC"
spe_clusters$PRECAST_cluster[spe_clusters$PRECAST_cluster == 5] <- "L6b"
spe_clusters$PRECAST_cluster[spe_clusters$PRECAST_cluster == 6] <- "L6a"
spe_clusters$PRECAST_cluster[spe_clusters$PRECAST_cluster == 4] <- "L5"
spe_clusters$PRECAST_cluster[spe_clusters$PRECAST_cluster == 2] <- "L3"
spe_clusters$PRECAST_cluster[spe_clusters$PRECAST_cluster == 1] <- "L2"
spe_clusters$PRECAST_cluster[spe_clusters$PRECAST_cluster == 9] <- "L1"

spe$PRECAST_cluster <- spe_clusters$PRECAST_cluster
dACC_layer <- "PRECAST_cluster"

make_plots <- function(factor, dACC_id){

    spe_1 <- spe[, which(spe$sample_id == dACC_id)]
    p1 <- make_escheR(spe_1) |>
        add_fill(var=factor, point_size = 3) |>
        add_ground(var=dACC_layer, stroke=0.5, point_size = 3) +
        scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", dACC_id)) +
        labs(color = "layer") +
        theme(
            legend.text = element_text(size = 5), # Smaller text size
            legend.title = element_text(size = 8), # Smaller title size
            legend.key.size = unit(0.5, "lines"), # Smaller key size
            legend.spacing.x = unit(0.3, 'cm'), # Adjust spacing between keys
            legend.spacing.y = unit(0.3, 'cm'))  # Adjust spacing between rows


    return(list(p1))
}

make_plots <- function(factor, dACC_id){

    spe_1 <- spe[, which(spe$sample_id == dACC_id)]
    p1 <- make_escheR(spe_1) |>
        add_fill(var=factor, point_size = 3) |>
        add_ground(var=dACC_layer, stroke=0.5, point_size = 3) +
        scale_color_manual(values = c(
            "L2" = "#E41A1C",   # Bright red
            "L3" = "#377EB8",   # Strong blue
            "L5" = "#4DAF4A",   # Vivid green
            "L6a" = "#984EA3",  # Purple
            "L6b" = "#FF7F00",  # Orange
            "WM1" = "#A65628",  # Rich brown
            "WM2" = "#FFFF99",  # Light yellow
            "WM-CC" = "#F781BF",# Pink
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


    return(list(p1))
}

#L2
plots_3 <- make_plots("NMF_3", "V12Y31-080_B1")

#L3
plots_10 <- make_plots("NMF_10", "V12N28-334_B1")

#L5
plots_26 <- make_plots("NMF_26", "V12Y31-080_C1")
plots_36 <- make_plots("NMF_36", "V12N28-331_C1")

#L6
plots_22 <- make_plots("NMF_22", "V12N28-334_C1")
plots_55 <- make_plots("NMF_55", "V12N28-332_C1")

#WM
plots_29 <- make_plots("NMF_29", "V12J03-002_A1")
plots_14 <- make_plots("NMF_14", "V12N28-331_D1")
plots_25 <- make_plots("NMF_25", "V12N28-332_D1")

pdf(here::here("plots", "13_NMF", "compare_NMF.pdf"), width = 10, height = 10)

grid.arrange(grobs = plots_3, nrow = 1)
grid.arrange(grobs = plots_10, nrow = 1)
grid.arrange(grobs = plots_26, nrow = 1)
grid.arrange(grobs = plots_36, nrow = 1)
grid.arrange(grobs = plots_22, nrow = 1)
grid.arrange(grobs = plots_55, nrow = 1)
grid.arrange(grobs = plots_29, nrow = 1)
grid.arrange(grobs = plots_14, nrow = 1)
grid.arrange(grobs = plots_25, nrow = 1)

dev.off()

