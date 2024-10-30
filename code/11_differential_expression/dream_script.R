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
library(BiocParallel)

param <- SnowParam(6, "SOCK", progressbar = TRUE)

load(
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

# create Layer_Tissue column
colData(sce_pseudo)$Layer_Tissue <- paste(colData(sce_pseudo)$Layer, colData(sce_pseudo)$Tissue, sep = "_")

combinations <- list(
    "L1_dACC_L1_DLPFC", "L1_dACC_L2_DLPFC", "L1_dACC_L3_DLPFC", "L1_dACC_L4_DLPFC", "L1_dACC_L5_DLPFC", "L1_dACC_L6_DLPFC", "L1_dACC_WM_DLPFC",
    "L2_dACC_L1_DLPFC", "L2_dACC_L2_DLPFC", "L2_dACC_L3_DLPFC", "L2_dACC_L4_DLPFC", "L2_dACC_L5_DLPFC", "L2_dACC_L6_DLPFC", "L2_dACC_WM_DLPFC",
    "L3_dACC_L1_DLPFC", "L3_dACC_L2_DLPFC", "L3_dACC_L3_DLPFC", "L3_dACC_L4_DLPFC", "L3_dACC_L5_DLPFC", "L3_dACC_L6_DLPFC", "L3_dACC_WM_DLPFC",
    "L5_dACC_L1_DLPFC", "L5_dACC_L2_DLPFC", "L5_dACC_L3_DLPFC", "L5_dACC_L4_DLPFC", "L5_dACC_L5_DLPFC", "L5_dACC_L6_DLPFC", "L5_dACC_WM_DLPFC",
    "L6a_dACC_L1_DLPFC", "L6a_dACC_L2_DLPFC", "L6a_dACC_L3_DLPFC", "L6a_dACC_L4_DLPFC", "L6a_dACC_L5_DLPFC", "L6a_dACC_L6_DLPFC", "L6a_dACC_WM_DLPFC",
    "L6b_dACC_L1_DLPFC", "L6b_dACC_L2_DLPFC", "L6b_dACC_L3_DLPFC", "L6b_dACC_L4_DLPFC", "L6b_dACC_L5_DLPFC", "L6b_dACC_L6_DLPFC", "L6b_dACC_WM_DLPFC",
    "WM_dACC_L1_DLPFC", "WM_dACC_L2_DLPFC", "WM_dACC_L3_DLPFC", "WM_dACC_L4_DLPFC", "WM_dACC_L5_DLPFC", "WM_dACC_L6_DLPFC", "WM_dACC_WM_DLPFC"
)

# Initialize a list to store results for all comparisons
results_list <- list()

# Loop through each combination to perform the DE analysis
for (comparison in combinations) {
    print(comparison)
    # Split the comparison into components
    parts <- unlist(strsplit(comparison, "_"))
    layer1 <- parts[1]  # e.g., L1
    tissue1 <- parts[2] # e.g., dACC
    layer2 <- parts[3]  # e.g., L1
    tissue2 <- parts[4] # e.g., DLPFC

    # Create the Layer_Tissue identifiers
    layer_tissue1 <- paste(layer1, tissue1, sep = "_")
    layer_tissue2 <- paste(layer2, tissue2, sep = "_")

    # Subset the data for the specific layers and tissues
    sce_pseudo_subset <- sce_pseudo[, colData(sce_pseudo)$Layer_Tissue %in% c(layer_tissue1, layer_tissue2)]

    # Proceed only if the subset has data
    if (ncol(sce_pseudo_subset) > 0) {

        dge <- DGEList(counts(sce_pseudo_subset))
        dge <- calcNormFactors(dge)

        # Update the design matrix with correct treatment of Layer_Tissue
        form <- "~ 0 + Layer_Tissue + (1 | Donor)"
        vobjDream <- voomWithDreamWeights(dge, form, colData(sce_pseudo_subset), BPPARAM=param)

        # Fit the model
        fit <- dream(vobjDream, form, colData(sce_pseudo_subset), BPPARAM=param)
        fit <- eBayes(fit)

        # find which coef in colnames(fit) starts with Layer_Tissue
        coef_name1 <- colnames(fit)[grep("Layer_Tissue", colnames(fit))]

        print(coef_name1)

        # Perform DE analysis
        results <- topTable(fit, coef = coef_name1, number = Inf)

        results_list[[comparison]] <- results
        # Count significant DE genes
        significant_genes <- sum(abs(results$z.std) > 1.645)
        print(paste("Significant DE genes for", comparison, ":", significant_genes))

    }
    else {
        print(paste("No data for", comparison))
    }
}

# Save the results
saveRDS(results_list, file = here("processed-data", "11_differential_expression", "dream_results.rds"))

# Define the layers for dACC and DLPFC
dACC_layers <- c("L1", "L2", "L3", "L5", "L6a", "L6b", "WM")
DLPFC_layers <- c("L1", "L2", "L3", "L4", "L5", "L6", "WM")

# Initialize a 7x7 matrix to store significant gene counts
heatmap_matrix <- matrix(0, nrow = length(dACC_layers), ncol = length(DLPFC_layers))
rownames(heatmap_matrix) <- dACC_layers
colnames(heatmap_matrix) <- DLPFC_layers

# Loop through each combination of dACC and DLPFC layers to count significant DE genes
for (dACC_layer in dACC_layers) {
    for (DLPFC_layer in DLPFC_layers) {
        comparison <- paste(dACC_layer, "dACC", DLPFC_layer, "DLPFC", sep = "_")

        # Extract results for the current comparison
        results <- results_list[[comparison]]

        # Count significant genes if results are available
        if (!is.null(results)) {
            heatmap_matrix[dACC_layer, DLPFC_layer] <- sum(results$adj.P.Val < 0.05)
        } else {
            heatmap_matrix[dACC_layer, DLPFC_layer] <- 0  # No data
        }
    }
}

# Convert the heatmap_matrix to a data frame for ggplot
library(reshape2)
heatmap_data <- melt(heatmap_matrix)

# Save the heatmap to a PDF file
pdf(here("plots", "11_differential_expression", "dream_heatmap.pdf"), width = 8, height = 6)

# Create the heatmap using ggplot2
library(ggplot2)
ggplot(heatmap_data, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "blue", high = "red") + # Adjust color gradient
    labs(title = "Number of Significant Genes by Layer Comparison",
         x = "DLPFC Layer",
         y = "dACC Layer") +
    theme_minimal() +                  # Use a minimal theme
    theme(axis.text.x = element_text(angle = 45, hjust = 1), # Adjust x-axis labels
          axis.text.y = element_text()) # Adjust y-axis labels

# Close the PDF device
dev.off()
