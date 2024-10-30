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
library(variancePartition)
library(edgeR)
library(reshape2)

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

# combined pseudobulked logcounts matrix from both tissues
spe_pseudo_dACC <-
    registration_pseudobulk(spe_dACC,
                            var_registration = "layer",
                            var_sample_id = "sample_id",
                            min_ncells = 10
    )
spe_pseudo_DLPFC <-
    registration_pseudobulk(spe_DLPFC,
                            var_registration = "BayesSpace_harmony_09",
                            var_sample_id = "sample_id",
                            min_ncells = 10
    )

spe_pseudo_dACC$Tissue <- "dACC"
spe_pseudo_DLPFC$Tissue <- "DLPFC"

intersected_genes <- intersect(rownames(spe_pseudo_dACC), rownames(spe_pseudo_DLPFC))
spe_pseudo_dACC <- spe_pseudo_dACC[intersected_genes, , drop = FALSE]
spe_pseudo_DLPFC <- spe_pseudo_DLPFC[intersected_genes, , drop = FALSE]

# save
save(
    spe_pseudo_dACC,
    spe_pseudo_DLPFC,
    file = here("processed-data", "11_differential_expression", "spe_pseudo_dACC_DLPFC.Rdata")
)


# create sample info for dACC with the counts matrix
info_dACC <- data.frame(
    Layer = spe_pseudo_dACC$layer,
    Donor = spe_pseudo_dACC$donor,
    Tissue = spe_pseudo_dACC$Tissue
)
rownames(info_dACC) <- colnames(counts(spe_pseudo_dACC))

# create sample info for DLPFC with the counts matrix
info_DLPFC <- data.frame(
    Layer = spe_pseudo_DLPFC$BayesSpace_harmony_09,
    Donor = spe_pseudo_DLPFC$donor,
    Tissue = spe_pseudo_DLPFC$Tissue
)
rownames(info_DLPFC) <- colnames(counts(spe_pseudo_DLPFC))

# combine gene expression matrix and sample info
countMatrix <- cbind(counts(spe_pseudo_dACC), counts(spe_pseudo_DLPFC))
info <- rbind(info_dACC, info_DLPFC)

# make spe
sce_pseudo <- SingleCellExperiment(
    assays = list(counts = countMatrix),
    colData = info
)

sce_pseudo <- logNormCounts(sce_pseudo)

pca <- prcomp(t(assays(sce_pseudo)$logcounts))
metadata(sce_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
pca_pseudo <- pca$x[, seq_len(20)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(sce_pseudo) <- list(PCA = pca_pseudo)

pdf(file = here("plots", "11_differential_expression",
                paste0("pseudobulk_PC_dACC_DLPFC.pdf")),
    width = 10, height = 10)

plotPCA(
    sce_pseudo,
    colour_by = "Layer",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(sce_pseudo)$PCA_var_explained
)

plotPCA(
    sce_pseudo,
    colour_by = "Donor",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(sce_pseudo)$PCA_var_explained
)

plotPCA(
    sce_pseudo,
    colour_by = "Tissue",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(sce_pseudo)$PCA_var_explained
)

vars <- getVarianceExplained(sce_pseudo,
                             variables = c("Layer","Donor", "Tissue")
)


plotExplanatoryVariables(vars)

dev.off()

################################## LIMMA ##################################

combinations <- list(
    "L1_dACC_L1_DLPFC", "L1_dACC_L2_DLPFC", "L1_dACC_L3_DLPFC", "L1_dACC_L4_DLPFC", "L1_dACC_L5_DLPFC", "L1_dACC_L6_DLPFC", "L1_dACC_WM_DLPFC",
    "L2_dACC_L1_DLPFC", "L2_dACC_L2_DLPFC", "L2_dACC_L3_DLPFC", "L2_dACC_L4_DLPFC", "L2_dACC_L5_DLPFC", "L2_dACC_L6_DLPFC", "L2_dACC_WM_DLPFC",
    "L3_dACC_L1_DLPFC", "L3_dACC_L2_DLPFC", "L3_dACC_L3_DLPFC", "L3_dACC_L4_DLPFC", "L3_dACC_L5_DLPFC", "L3_dACC_L6_DLPFC", "L3_dACC_WM_DLPFC",
    "L5_dACC_L1_DLPFC", "L5_dACC_L2_DLPFC", "L5_dACC_L3_DLPFC", "L5_dACC_L4_DLPFC", "L5_dACC_L5_DLPFC", "L5_dACC_L6_DLPFC", "L5_dACC_WM_DLPFC",
    "L6a_dACC_L1_DLPFC", "L6a_dACC_L2_DLPFC", "L6a_dACC_L3_DLPFC", "L6a_dACC_L4_DLPFC", "L6a_dACC_L5_DLPFC", "L6a_dACC_L6_DLPFC", "L6a_dACC_WM_DLPFC",
    "L6b_dACC_L1_DLPFC", "L6b_dACC_L2_DLPFC", "L6b_dACC_L3_DLPFC", "L6b_dACC_L4_DLPFC", "L6b_dACC_L5_DLPFC", "L6b_dACC_L6_DLPFC", "L6b_dACC_WM_DLPFC",
    "WM_dACC_L1_DLPFC", "WM_dACC_L2_DLPFC", "WM_dACC_L3_DLPFC", "WM_dACC_L4_DLPFC", "WM_dACC_L5_DLPFC", "WM_dACC_L6_DLPFC", "WM_dACC_WM_DLPFC"
)

results_list <- list()

# create Layer_Tissue column
colData(sce_pseudo)$Layer_Tissue <- paste(colData(sce_pseudo)$Layer, colData(sce_pseudo)$Tissue, sep = "_")

for (comparison in combinations) {
    print(comparison)
    parts <- unlist(strsplit(comparison, "_"))
    layer1 <- parts[1]
    tissue1 <- parts[2]
    layer2 <- parts[3]
    tissue2 <- parts[4]

    layer_tissue1 <- paste(layer1, tissue1, sep = "_")
    layer_tissue2 <- paste(layer2, tissue2, sep = "_")

    sce_pseudo_subset <- sce_pseudo[, colData(sce_pseudo)$Layer_Tissue %in% c(layer_tissue1, layer_tissue2)]

    if (ncol(sce_pseudo_subset) > 0) {
        mat_formula <- "~ Donor + Layer_Tissue"
        mod <- model.matrix(as.formula(mat_formula), data = colData(sce_pseudo_subset))

        fit <- lmFit(logcounts(sce_pseudo_subset), mod)
        fit <- eBayes(fit)

        # find which coef in colnames(fit) starts with Layer_Tissue
        coef_name1 <- colnames(fit)[grep("Layer_Tissue", colnames(fit))]

        print(coef_name1)

        results <- topTable(fit, coef = coef_name1, adjust = "fdr", number = Inf)
        results_list[[comparison]] <- results
        significant_genes <- sum(results$adj.P.Val < 0.05)
        print(paste("Significant DE genes for", comparison, ":", significant_genes))


    }
    else {
        print(paste("No data for", comparison))
    }
}

dACC_layers <- c("L1", "L2", "L3", "L5", "L6a", "L6b", "WM")
DLPFC_layers <- c("L1", "L2", "L3", "L4", "L5", "L6", "WM")

heatmap_matrix <- matrix(0, nrow = length(dACC_layers), ncol = length(DLPFC_layers))
rownames(heatmap_matrix) <- dACC_layers
colnames(heatmap_matrix) <- DLPFC_layers

for (dACC_layer in dACC_layers) {
    for (DLPFC_layer in DLPFC_layers) {
        comparison <- paste(dACC_layer, "dACC", DLPFC_layer, "DLPFC", sep = "_")

        results <- results_list[[comparison]]

        if (!is.null(results)) {
            heatmap_matrix[dACC_layer, DLPFC_layer] <- sum(results$adj.P.Val < 0.05)
        } else {
            heatmap_matrix[dACC_layer, DLPFC_layer] <- 0  # No data
        }
    }
}

heatmap_data <- melt(heatmap_matrix)

pdf(here("plots", "11_differential_expression", "limma_heatmap.pdf"), width = 8, height = 6)

ggplot(heatmap_data, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "blue", high = "red") +
    labs(title = "Number of Significant Genes by Layer Comparison",
         x = "DLPFC Layer",
         y = "dACC Layer") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text())

dev.off()
