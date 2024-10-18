setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library(dplyr)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(EnhancedVolcano)
library(scater)
library(here)
library(spatialLIBD)
library(DeconvoBuddies)

# load spe for k=9 without WM-CC
load(
    file = here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata")
)
spe_dACC <- spe
spe_dACC$donor <- spe_dACC$brnum

# load DLPFC spe
spe_DLPFC <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium")
spe_DLPFC$donor <- spe_DLPFC$subject

# create spatial labels for DLPFC_30
spe_DLPFC$BayesSpace_harmony_09[spe_DLPFC$BayesSpace_harmony_09 == 3] <- "L2"
spe_DLPFC$BayesSpace_harmony_09[spe_DLPFC$BayesSpace_harmony_09 == 8] <- "L4"
spe_DLPFC$BayesSpace_harmony_09[spe_DLPFC$BayesSpace_harmony_09 == 7] <- "L6"
spe_DLPFC$BayesSpace_harmony_09[spe_DLPFC$BayesSpace_harmony_09 == 5] <- "L3"
spe_DLPFC$BayesSpace_harmony_09[spe_DLPFC$BayesSpace_harmony_09 == 6] <- "WM"
spe_DLPFC$BayesSpace_harmony_09[spe_DLPFC$BayesSpace_harmony_09 == 4] <- "L5"
spe_DLPFC$BayesSpace_harmony_09[spe_DLPFC$BayesSpace_harmony_09 == 2] <- "L1"
spe_DLPFC$BayesSpace_harmony_09[spe_DLPFC$BayesSpace_harmony_09 == 1] <- "meninges"
spe_DLPFC$BayesSpace_harmony_09[spe_DLPFC$BayesSpace_harmony_09 == 9] <- "WM"
# remove meninges
spe_DLPFC <- spe_DLPFC[ , which(spe_DLPFC$BayesSpace_harmony_09 != "meninges")]

brains <- unique(colData(spe_dACC)$donor)

intersected_genes <- intersect(rownames(spe_dACC), rownames(spe_DLPFC))
spe_dACC <- spe_dACC[intersected_genes, , drop = FALSE]
spe_DLPFC <- spe_DLPFC[intersected_genes, , drop = FALSE]

all_corr_list <- list()
all_high_corr_list <- list()

for (brain in brains) {
    print(paste("Processing donor:", brain))

    spe_dACC_Br <- spe_dACC[, colData(spe_dACC)$donor == brain]
    spe_DLPFC_Br <- spe_DLPFC[, colData(spe_DLPFC)$donor == brain]

    markers_dACC_Br <- scran::findMarkers(
        spe_dACC_Br,
        groups = colData(spe_dACC_Br)$layer,
        assay.type = "logcounts",
        sorted = T,
        test.type = "t"
    )

    markers_DLPFC_Br <- scran::findMarkers(
        spe_DLPFC_Br,
        groups = colData(spe_DLPFC_Br)$BayesSpace_harmony_09,
        assay.type = "logcounts",
        sorted = T,
        test.type = "t"
    )

    corr_list <- list()
    high_corr_list <- list()

    for (layer in names(markers_DLPFC_Br)) {
        print(layer)

        # variable to keep track of the highest correlation for the current layer
        max_corr <- -Inf
        best_dACC_layer <- NULL

        # top 100 marker genes for DLPFC layer
        markers_DLPFC_Br[[layer]] <- markers_DLPFC_Br[[layer]][markers_DLPFC_Br[[layer]]$Top <= 100, ]

        for (dACC_layer in names(markers_dACC_Br)) {
            print(dACC_layer)

            # intersect genes from both datasets, DLPFC dataset only has 100
            i <- intersect(rownames(markers_dACC_Br[[dACC_layer]]), rownames(markers_DLPFC_Br[[layer]]))
            markers_dACC_Br[[dACC_layer]] <- markers_dACC_Br[[dACC_layer]][i, ]

            # match indices between the two datasets
            markers_dACC_Br[[dACC_layer]] <- markers_dACC_Br[[dACC_layer]][rownames(markers_DLPFC_Br[[layer]]), ]
            corr <- cor(markers_dACC_Br[[dACC_layer]]$summary.logFC, markers_DLPFC_Br[[layer]]$summary.logFC, use = "complete.obs")
            corr_list[[paste(layer, dACC_layer, sep = "_")]] <- corr
            print(corr)

            if (corr > max_corr) {
                max_corr <- corr
                best_dACC_layer <- dACC_layer
                }
        }

        high_corr_list[[layer]] <- list(best_dACC_layer = best_dACC_layer, max_corr = max_corr)
    }

    all_corr_list[[brain]] <- corr_list
    all_high_corr_list[[brain]] <- high_corr_list
}

print(all_high_corr_list)

results_table <- list()

for (brain in names(all_high_corr_list)) {
    high_corr_list <- all_high_corr_list[[brain]]

    for (layer in names(high_corr_list)) {
        best_dACC_layer <- high_corr_list[[layer]]$best_dACC_layer
        max_corr <- high_corr_list[[layer]]$max_corr

        results_table[[length(results_table) + 1]] <- data.frame(
            Donor = brain,
            Layer = layer,
            Best_dACC_Layer = best_dACC_layer,
            Max_Correlation = max_corr
        )
    }
}

results_table_df <- do.call(rbind, results_table)

pdf(here("plots", "11_differential_expression", "donor_DLPFC_dACC_correlation_boxplot.pdf"))
ggplot(results_table_df, aes(x = Layer, y = Max_Correlation, fill = Layer)) +
    geom_boxplot() +
    labs(title = "Box Plot of Max Correlation by DLPFC Layer",
         x = "DLPFC Layer",
         y = "Max Correlation") +
    theme_minimal() +
    theme(legend.position = "none")
dev.off()

pdf(here("plots", "11_differential_expression", "donor_DLPFC_dACC_correlation_boxplot_batch.pdf"))
ggplot(results_table_df, aes(x = Donor, y = Max_Correlation, fill = Donor)) +
    geom_boxplot() +
    labs(title = "Box Plot of Correlations by Donor",
         x = "Donor",
         y = "Max Correlation") +
    theme_minimal() +
    theme(legend.title = element_blank())
dev.off()
