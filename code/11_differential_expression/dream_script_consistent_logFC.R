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
library(reshape2)

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

results_list <- list()

for (comparison in combinations) {
    print(comparison)
    parts <- unlist(strsplit(comparison, "_"))
    layer1 <- parts[1]
    tissue1 <- parts[2]
    layer2 <- parts[3]
    tissue2 <- parts[4]

    layer_tissue1 <- paste(layer1, tissue1, sep = "_")
    layer_tissue2 <- paste(layer2, tissue2, sep = "_")

    # Subset the data
    sce_pseudo_subset <- sce_pseudo[, colData(sce_pseudo)$Layer_Tissue %in% c(layer_tissue1, layer_tissue2)]

    if (ncol(sce_pseudo_subset) > 0) {
        # Explicitly set the factor levels to ensure dACC is the reference
        colData(sce_pseudo_subset)$Layer_Tissue <- factor(
            colData(sce_pseudo_subset)$Layer_Tissue,
            levels = c(layer_tissue1, layer_tissue2)  # Ensures dACC is first
        )

        # Perform the analysis
        dge <- DGEList(counts(sce_pseudo_subset))
        dge <- calcNormFactors(dge)

        form <- "~ Layer_Tissue + (1 | Donor)"
        vobjDream <- voomWithDreamWeights(dge, form, colData(sce_pseudo_subset), BPPARAM=param)

        fit <- dream(vobjDream, form, colData(sce_pseudo_subset), BPPARAM=param)
        fit <- eBayes(fit)

        results <- topTable(fit, number = Inf)

        # Ensure logFC is consistent
        # Positive logFC means upregulated in dACC
        if (levels(colData(sce_pseudo_subset)$Layer_Tissue)[1] != layer_tissue1) {
            stop("Factor levels were not set correctly")
        }

        results_list[[comparison]] <- results
        significant_genes <- sum(abs(results$z.std) > 1.645)
        print(paste("Significant DE genes for", comparison, ":", significant_genes))
    } else {
        print(paste("No data for", comparison))
    }
}


# save the results
saveRDS(results_list, file = here("processed-data", "11_differential_expression", "dream_results_consistent_logFC.rds"))
