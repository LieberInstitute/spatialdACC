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

    sce_pseudo_subset <- sce_pseudo[, colData(sce_pseudo)$Layer_Tissue %in% c(layer_tissue1, layer_tissue2)]

    if (ncol(sce_pseudo_subset) > 0) {
        # Explicitly set the factor levels to ensure dACC is the reference
        colData(sce_pseudo_subset)$Layer_Tissue <- factor(
            colData(sce_pseudo_subset)$Layer_Tissue,
            levels = c(layer_tissue1, layer_tissue2)
        )

        dge <- DGEList(counts(sce_pseudo_subset))
        dge <- calcNormFactors(dge)

        form <- "~ Layer_Tissue + (1 | Donor)"
        vobjDream <- voomWithDreamWeights(dge, form, colData(sce_pseudo_subset), BPPARAM=param)

        fit <- dream(vobjDream, form, colData(sce_pseudo_subset), BPPARAM=param)
        fit <- eBayes(fit)

        results <- topTable(fit, number = Inf)

        # positive logFC means upregulated in DLPFC
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
            df <- results %>%
                filter(abs(z.std) > 1.645 & abs(logFC) > log2(1.5))
            heatmap_matrix[dACC_layer, DLPFC_layer] <- dim(df)[1]
        } else {
            heatmap_matrix[dACC_layer, DLPFC_layer] <- 0
        }
    }
}

heatmap_data <- melt(heatmap_matrix)

# remove rows that have "WM"
heatmap_data <- heatmap_data[!grepl("WM", heatmap_data$Var1),]
heatmap_data <- heatmap_data[!grepl("WM", heatmap_data$Var2),]

pdf(here("plots", "11_differential_expression", "dream_heatmap_noWM.pdf"), width = 8, height = 6)

ggplot(heatmap_data, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    geom_text(aes(label=value), color="black", size=2) +
    scale_fill_gradient(low = "grey", high = "red") +
    labs(title = "Number of Significant Genes by Layer Comparison",
         x = "DLPFC Layer",
         y = "dACC Layer",
         caption = "using logFC threshold of 1.5 with dream") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text())

dev.off()


# add up total number of significant genes for each Var2 (DLPFC layer)
heatmap_data %>%
    group_by(Var2) %>%
    summarise(total_significant_genes = sum(value)) %>%
    arrange(desc(total_significant_genes))

#1 L1                      25102
#2 L4                      24593
#3 WM                      23989
#4 L5                      22753
#5 L3                      21902
#6 L2                      21516
#7 L6                      20710

# standardize the heatmap data by dividing by the total number of significant genes for each DLPFC layer
heatmap_data_DLPFC <- heatmap_data %>%
    group_by(Var2) %>%
    mutate(value = value / sum(value))

# round to two decimal places
heatmap_data_DLPFC$value <- round(heatmap_data_DLPFC$value, 2)

pdf(here("plots", "11_differential_expression", "dream_heatmap_DLPFC_stand.pdf"), width = 8, height = 6)

ggplot(heatmap_data_DLPFC, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    geom_text(aes(label=value), color="black", size=2) +
    scale_fill_gradient(low = "grey", high = "red") +
    labs(title = "Number of Significant Genes by Layer Comparison",
         x = "DLPFC Layer",
         y = "dACC Layer",
         caption = "standardized by total DLPFC sig genes, removing WM") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text())

dev.off()

heatmap_data %>%
    group_by(Var1) %>%
    summarise(total_significant_genes = sum(value)) %>%
    arrange(desc(total_significant_genes))
#1 L1                      21562
#2 L6b                     21019
#3 L2                      13583
#4 L5                      11917
#5 L6a                     11561
#6 L3                      11472

# standardize the heatmap data by dividing by the total number of significant genes for each dACC layer
heatmap_data_dACC <- heatmap_data %>%
    group_by(Var1) %>%
    mutate(value = value / sum(value))

# round to two decimal places
heatmap_data_dACC$value <- round(heatmap_data_dACC$value, 2)

pdf(here("plots", "11_differential_expression", "dream_heatmap_dACC_stand.pdf"), width = 8, height = 6)

ggplot(heatmap_data_dACC, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    geom_text(aes(label=value), color="black", size=2) +
    scale_fill_gradient(low = "blue", high = "white") +
    labs(title = "Number of Significant Genes by Layer Comparison",
         x = "DLPFC Layer",
         y = "dACC Layer",
         caption = "standardized by total dACC sig genes, removing WM") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text())

dev.off()

# make volcano plots for all comparisons
for (comparison in names(results_list)) {
    results <- results_list[[comparison]]
    print(comparison)

    results$gene_name <- rowData(spe_pseudo_dACC)[rownames(results),]$gene_name

    p <- ggplot(results, aes(x = logFC, y = -log10(P.Value))) +
        geom_point(aes(color = ifelse(abs(z.std) > 1.645, "red", "black")), alpha = 0.5) +
        scale_color_identity() +
        geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        labs(title = paste("Volcano Plot for", comparison),
             x = "log2 Fold Change",
             y = "-log10 P-value",
             caption = "using dream, positive x-axis means upregulated in DLPFC") +
        theme_bw() +
        geom_text_repel(data = results %>% filter(abs(z.std) > 1.645 & abs(logFC) > log2(1.5)),
                        aes(label = gene_name),
                        box.padding = 0.5,
                        point.padding = 0.5,
                        segment.color = "grey50",
                        segment.size = 0.5,
                        segment.alpha = 0.5,
                        size = 2)

    Sys.sleep(1)

    pdf(here("plots", "11_differential_expression", "dream_volcano",paste0("dream_volcano_", comparison, ".pdf")), width = 8, height = 6)
    print(p)
    dev.off()
}


