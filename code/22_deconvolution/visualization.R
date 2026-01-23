setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library(here)
library(SpatialExperiment)
library(spacexr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(escheR)
library(ComplexHeatmap)
library(circlize)

# set directories
plots_dir <- here("plots", "22_deconvolution")
processed_dir <- here("processed-data", "22_deconvolution")


#set up for dACC into dACC
##################################################################################################################################
# set file prefix, sample id var, layer colors
file_prefix <- "dACC_dACC"
sample_id_var <- "sample_id"
layer_var <- "layer"
layer_order <- c("L1", "L2", "L3", "L5", "L6a", "L6b", "WM")

layer_colors <- c(
    "L1" = "#F0027F",
    "L2" = "#377EB8",
    "L3" = "#4DAF4A",
    "L5" = "#FFD700",
    "L6a" = "#FFC18A",
    "L6b" = "#c46200",
    "WM" = "#1A1A1A"
)

# load desired spe
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata"))
spe <- spe[!duplicated(rownames(spe)), ]

# load corresponding res file
res <- readRDS(here(processed_dir, paste0("rctd_results_",file_prefix,".rds")))

##################################################################################################################################

#set up for dACC into dlPFC
##################################################################################################################################
# set file prefix, sample id var, layer colors
file_prefix <- "dACC_dlPFC"
sample_id_var <- "sample_id"
layer_var <- "BayesSpace_harmony_09"
layer_order <- c("L1", "L2", "L3", "L4", "L5", "L6", "WM")

layer_colors <- c(
    "L2" = "#377EB8",
    "L3" = "#4DAF4A",
    "L5" = "#FFD700",
    "L4" = "#984EA3",
    "L6" = "#FF7F00",
    "WM" = "#1A1A1A",
    "L1" = "#F0027F"
)

# load desired spe
spe <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium")
spe_DLPFC_30.temp <- spe
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 3] <- "L2"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 8] <- "L4"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 7] <- "L6"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 5] <- "L3"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 6] <- "WM"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 4] <- "L5"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 2] <- "L1"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 1] <- "meninges"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 9] <- "WM"
spe_DLPFC_30.temp <- spe_DLPFC_30.temp[ , which(spe_DLPFC_30.temp$BayesSpace_harmony_09 != "meninges")]
spe <- spe_DLPFC_30.temp
spe <- spe[!duplicated(rownames(spe)), ]
colnames(spe) <- paste0(colnames(spe),colData(spe)$sample_id)

# load corresponding res file
res <- readRDS(here(processed_dir, paste0("rctd_results_",file_prefix,".rds")))

##################################################################################################################################


#set up for dlPFC into dlPFC
##################################################################################################################################
# set file prefix, sample id var, layer colors
file_prefix <- "dlPFC_dlPFC"
sample_id_var <- "sample_id"
layer_var <- "BayesSpace_harmony_09"
layer_order <- c("L1", "L2", "L3", "L4", "L5", "L6", "WM")

layer_colors <- c(
    "L2" = "#377EB8",
    "L3" = "#4DAF4A",
    "L5" = "#FFD700",
    "L4" = "#984EA3",
    "L6" = "#FF7F00",
    "WM" = "#1A1A1A",
    "L1" = "#F0027F"
)

# load desired spe
spe <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium")
spe_DLPFC_30.temp <- spe
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 3] <- "L2"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 8] <- "L4"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 7] <- "L6"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 5] <- "L3"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 6] <- "WM"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 4] <- "L5"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 2] <- "L1"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 1] <- "meninges"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 9] <- "WM"
spe_DLPFC_30.temp <- spe_DLPFC_30.temp[ , which(spe_DLPFC_30.temp$BayesSpace_harmony_09 != "meninges")]
spe <- spe_DLPFC_30.temp
spe <- spe[!duplicated(rownames(spe)), ]
colnames(spe) <- paste0(colnames(spe),colData(spe)$sample_id)

# load corresponding res file
res <- readRDS(here(processed_dir, paste0("rctd_results_",file_prefix,".rds")))

##################################################################################################################################

#set up for dlPFC into dACC
##################################################################################################################################
# set file prefix, sample id var, layer colors
file_prefix <- "dlPFC_dACC"
sample_id_var <- "sample_id"
layer_var <- "layer"
layer_order <- c("L1", "L2", "L3", "L5", "L6a", "L6b", "WM")

layer_colors <- c(
    "L1" = "#F0027F",
    "L2" = "#377EB8",
    "L3" = "#4DAF4A",
    "L5" = "#FFD700",
    "L6a" = "#FFC18A",
    "L6b" = "#c46200",
    "WM" = "#1A1A1A"
)

# load desired spe
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata"))
spe <- spe[!duplicated(rownames(spe)), ]
colnames(spe) <- paste0(colnames(spe),colData(spe)$sample_id)

# load corresponding res file
res <- readRDS(here(processed_dir, paste0("rctd_results_",file_prefix,".rds")))

##################################################################################################################################


# Extract weights matrix and transpose (so spots are rows, cell types are columns)
weights_mat <- t(assays(res)$weights)

# Create weights dataframe
weights <- as.data.frame(as.matrix(weights_mat))

# Match spots between res and spe
common_spots <- intersect(rownames(weights), colnames(spe))
length(common_spots)

# Subset to common spots
weights <- weights[common_spots, ]
spe_sub <- spe[, common_spots]

# Add spatial coordinates and sample info
coords <- spatialCoords(spe_sub)
weights$x <- coords[, 1]
weights$y <- coords[, 2]
weights$sample_id <- colData(spe_sub)[,sample_id_var]

# Get cell types and samples
cell_types <- colnames(weights_mat)
samples <- unique(weights$sample_id)

pdf(file.path(plots_dir, paste0(file_prefix,"_rctd_celltype_weights_by_sample.pdf")), width = 16, height = 12)

for (samp in samples) {
    weights_sub <- weights[weights$sample_id == samp, ]

    weights_long <- weights_sub %>%
        pivot_longer(
            cols = all_of(cell_types),
            names_to = "cell_type",
            values_to = "weight"
        )

    p <- ggplot(weights_long, aes(x = x, y = y, color = weight)) +
        geom_point(size = 0.5) +
        scale_color_viridis_c() +
        scale_y_reverse() +
        facet_wrap(~cell_type, ncol = 5) +
        labs(title = paste0("Sample: ", samp)) +
        theme_minimal() +
        coord_fixed()

    print(p)
}

dev.off()

# help visualize for dlPFC -> dlPFC because the cell type signal is so weak
# Calculate max weight for each cell type
max_weights <- apply(weights[, cell_types], 2, max, na.rm = TRUE)

pdf(file.path(plots_dir, paste0(file_prefix, "_rctd_celltype_weights_by_celltype.pdf")), width = 20, height = 16)

for (ct in cell_types) {
    # Get max weight for this cell type
    max_wt <- max_weights[ct]

    p <- ggplot(weights, aes(x = x, y = y, color = .data[[ct]])) +
        geom_point(size = 0.3) +
        scale_color_viridis_c(limits = c(0, max_wt), name = "Weight") +
        scale_y_reverse() +
        facet_wrap(~sample_id) +
        labs(title = paste0("Cell Type: ", ct)) +
        theme_minimal() +
        theme(aspect.ratio = 1)  # Use aspect.ratio instead of coord_fixed()

    print(p)
}

dev.off()

# Add weights to colData of spe_sub
for (ct in colnames(weights_mat)) {
    colData(spe_sub)[[ct]] <- weights[[ct]]
}

pdf(file.path(plots_dir, paste0(file_prefix,"_rctd_celltype_weights_escheR.pdf")), width = 16, height = 12)

for (samp in samples) {
    idx <- colData(spe_sub)$sample_id == samp
    spe_samp <- spe_sub[, idx]

    plot_list <- list()

    for (ct in cell_types) {
        p <- make_escheR(spe_samp) |>
            add_ground(var = layer_var, point_size = 0.8) |>
            add_fill(var = ct, stroke = 0, point_size = 0.6) +  # stroke = 0 removes outline
            scale_color_manual(values = layer_colors) +
            scale_fill_gradient(low = "white", high = "#54278F") +
            ggtitle(ct) +
            theme(
                legend.position = "none",
                plot.title = element_text(size = 10, hjust = 0.5)
            )

        plot_list[[ct]] <- p
    }

    combined <- patchwork::wrap_plots(plot_list, ncol = 5) +
        patchwork::plot_annotation(title = paste0("Sample: ", samp))

    print(combined)
}

dev.off()

# Extract weights matrix and transpose
weights_mat <- t(assays(res)$weights)
weights <- as.data.frame(as.matrix(weights_mat))

# Match spots
common_spots <- intersect(rownames(weights), colnames(spe))
weights <- weights[common_spots, ]
spe_sub <- spe[, common_spots]

# Add layer information
weights$layer <- colData(spe_sub)[, layer_var]

# Calculate mean weight per cell type per layer
mean_weights <- weights %>%
    group_by(layer) %>%
    summarise(across(all_of(colnames(weights_mat)), mean, na.rm = TRUE)) %>%
    as.data.frame()

#remove Excit_L5 for dlPFC_dACC - all zeroes
#mean_weights$Excit_L5 <- NULL

rownames(mean_weights) <- mean_weights$layer
mean_weights$layer <- NULL

heatmap_mat <- as.matrix(mean_weights)
heatmap_mat <- t(heatmap_mat)

# Order layers anatomically
layer_order <- layer_order[layer_order %in% colnames(heatmap_mat)]
heatmap_mat <- heatmap_mat[, layer_order]

# Scale by row
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

# Order cell types by their peak layer (to create diagonal)
# Find which layer has the max value for each cell type
peak_layer <- apply(heatmap_mat_scaled, 1, which.max)

# Create ordering based on peak layer, then by the value at that peak (descending)
cell_type_order <- names(sort(peak_layer +
                                  apply(heatmap_mat_scaled, 1, max) / (max(heatmap_mat_scaled) + 1)))

# Reorder matrix
heatmap_mat_scaled <- heatmap_mat_scaled[cell_type_order, ]

# Column annotation
col_anno <- HeatmapAnnotation(
    Layer = colnames(heatmap_mat_scaled),
    col = list(Layer = layer_colors[colnames(heatmap_mat_scaled)]),
    show_legend = FALSE
)

# Plot
pdf(file.path(plots_dir, paste0(file_prefix, "_rctd_celltype_layer_heatmap.pdf")), width = 8, height = 10)

ht <- Heatmap(
    heatmap_mat_scaled,
    name = "Scaled\nWeight",
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    top_annotation = col_anno,
    cluster_columns = FALSE,
    cluster_rows = FALSE,  # No clustering, use diagonal order
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_side = "left",
    column_names_rot = 45,
    column_title = "Spatial Domain",
    row_title = "Cell Type"
)

draw(ht)

dev.off()
