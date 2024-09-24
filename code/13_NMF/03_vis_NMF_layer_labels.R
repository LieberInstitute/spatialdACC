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
load(file = here("processed-data", "13_NMF", "spe_NMF.Rdata"))

dACC_layer <- "layer"

make_plots <- function(factor, dACC_id){

    spe_1 <- spe[, which(spe$sample_id == dACC_id)]
    p1 <- make_escheR(spe_1) |>
        add_fill(var=factor, point_size = 3) |>
        add_ground(var=dACC_layer, stroke=0.5, point_size = 3) +
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
        labs(color = "layer", caption = "single nucleus dACC NMF") +
        theme(
            legend.text = element_text(size = 5), # Smaller text size
            legend.title = element_text(size = 8), # Smaller title size
            legend.key.size = unit(0.5, "lines"), # Smaller key size
            legend.spacing.x = unit(0.3, 'cm'), # Adjust spacing between keys
            legend.spacing.y = unit(0.3, 'cm'))  # Adjust spacing between rows


    return(list(p1))
}

#L1
plots_21 <- make_plots("nmf21", "V12N28-331_C1")
plots_64 <- make_plots("nmf64", "V12N28-332_B1")
plots_65 <- make_plots("nmf65", "V12N28-334_C1")

#L2
plots_3 <- make_plots("nmf3", "V12Y31-080_B1")

#L3
plots_11 <- make_plots("nmf11", "V12N28-332_D1")

#L5a
plots_38 <- make_plots("nmf38", "V12Y31-080_C1")

#L5b
plots_61 <- make_plots("nmf61", "V12N28-332_D1")

#L6a
plots_15 <- make_plots("nmf15", "V12N28-332_A1")

#kind of L6b and a
plots_35 <- make_plots("nmf35", "V12N28-331_A1")

#WM
plots_6 <- make_plots("nmf6", "V12J03-002_A1")
plots_13 <- make_plots("nmf13", "V12N28-331_D1")
plots_27 <- make_plots("nmf27", "V12N28-332_D1")
plots_33 <- make_plots("nmf33", "V12N28-334_D1")

#misc
plots_5 <- make_plots("nmf5", "V12N28-332_A1")
plots_18 <- make_plots("nmf18", "V12N28-331_B1")
plots_49 <- make_plots("nmf49", "V12J03-002_C1")
plots_59 <- make_plots("nmf59", "V12N28-334_B1")

pdf(here::here("plots", "13_NMF", "compare_NMF.pdf"), width = 10, height = 10)

grid.arrange(grobs = plots_65, nrow = 1)
grid.arrange(grobs = plots_3, nrow = 1)
grid.arrange(grobs = plots_11, nrow = 1)
grid.arrange(grobs = plots_38, nrow = 1)
grid.arrange(grobs = plots_61, nrow = 1)
grid.arrange(grobs = plots_15, nrow = 1)
grid.arrange(grobs = plots_35, nrow = 1)
grid.arrange(grobs = plots_6, nrow = 1)
grid.arrange(grobs = plots_13, nrow = 1)
grid.arrange(grobs = plots_27, nrow = 1)
grid.arrange(grobs = plots_33, nrow = 1)
grid.arrange(grobs = plots_5, nrow = 1)
grid.arrange(grobs = plots_18, nrow = 1)
grid.arrange(grobs = plots_49, nrow = 1)
grid.arrange(grobs = plots_59, nrow = 1)

dev.off()

