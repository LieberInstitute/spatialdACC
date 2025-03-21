setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("scater")
    library("spatialLIBD")
    library("dplyr")
    library('EnhancedVolcano')
    library('ggrepel')
})

nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", 9)
load(file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_DE", paste0(nnSVG_precast_name,".Rdata"))
)
modeling_results_dACC <- modeling_results

load(file = here("processed-data", "12_spatial_registration",paste0("DLPFC_30_DE",".Rdata")))
modeling_results_DLPFC <- spe_modeling_results

# intersect genes
genes_dACC <- rownames(modeling_results_dACC[["enrichment"]])
genes_DLPFC <- rownames(modeling_results_DLPFC[["enrichment"]])
genes_intersect <- intersect(genes_dACC, genes_DLPFC)
dACC_results <- modeling_results_dACC[["enrichment"]][genes_intersect,]
DLPFC_results <- modeling_results_DLPFC[["enrichment"]][genes_intersect,]

pdf(file = here("plots", "11_differential_expression", "layer_markers_dACC_DLPFC.pdf"), width = 15, height = 15)

# repeat this in a loop for layers L1, L2, L3, L5
layers <- c("L1", "L2", "L3", "L5")
p_list <- list()

for (layer in layers) {
    genes <- dACC_results[, "gene"]
    dACC_logFC <- dACC_results[, paste0("logFC_", layer)]
    DLPFC_logFC <- DLPFC_results[, paste0("logFC_", layer)]
    corr_layer <- cor.test(dACC_logFC, DLPFC_logFC)

    # create data frame
    df_layer <- data.frame(
        gene = genes,
        dACC_logFC = dACC_logFC,
        DLPFC_logFC = DLPFC_logFC
    )

    df_layer$color <- with(df_layer, ifelse(
        dACC_logFC > 1.25 & DLPFC_logFC > 1.25, "DLPFC & dACC enriched",
        ifelse(dACC_logFC > 1.25 & DLPFC_logFC < 1.25, "dACC enriched", "Neither")
    ))

    # write to csv file
    # genes that are either "DLPFC & dACC enriched" or "dACC enriched"
    #write.csv(df_layer[df_layer$color != "Neither",],
    #          here("processed-data", "11_differential_expression", "novel_markers", paste0("dACC_DLPFC_", layer, ".csv")))

    if(layer == "L1"){
        genes_to_label <- c("RELN", "MSX1", "VIM","HBB","NTS","HBA1")
    }
    if(layer == "L2"){
        genes_to_label <- c("STXBP6", "LAMP5", "KCTD4", "ARHGAP4", "RSPO2", "CCNO","C1QL2")
    }
    if(layer == "L3"){
        genes_to_label <- c("LINC01007", "ADCYAP1")
    }
    if(layer == "L5"){
        genes_to_label <- c("PCP4", "TRABD2A", "MEPE", "CD24", "CD52", "FDPS", "DRD5", "GYG1", "ITGB1BP1")
    }

    # plot
    p <- ggplot(df_layer, aes(x = dACC_logFC, y = DLPFC_logFC)) +
        geom_point(aes(color = color), size=0.7) +
        scale_color_manual(values = c("DLPFC & dACC enriched" = "blue", "dACC enriched" = "red", "Neither" = "grey")) +
        geom_text_repel(aes(label = ifelse(gene %in% genes_to_label, gene, "")), size = 3, max.overlaps = Inf) +
        geom_smooth(method = "lm", se = FALSE, color = "blue") +
        geom_vline(xintercept = 1.25, color = "black", linetype = "dashed") +
        geom_hline(yintercept = 1.25, color = "black", linetype = "dashed") +
        geom_vline(xintercept = -1.25, color = "black", linetype = "dashed") +
        geom_hline(yintercept = -1.25, color = "black", linetype = "dashed") +
        xlim(-5, 5) +
        ylim(-5, 5) +
        labs(
            title = paste0("Layer: ", layer),
            x = "dACC logFC",
            y = "DLPFC logFC",
            subtitle = paste0("Pearson's r: ", round(corr_layer$estimate, 2))
        )

    print(p)
    p_list[[layer]] <- p
}

layers <- c("L6a", "L6b")
for (layer in layers) {
    genes <- dACC_results[, "gene"]
    dACC_logFC <- dACC_results[, paste0("logFC_", layer)]
    DLPFC_logFC <- DLPFC_results[, paste0("logFC_", "L6")]
    corr_layer <- cor.test(dACC_logFC, DLPFC_logFC)

    # create data frame
    df_layer <- data.frame(
        gene = genes,
        dACC_logFC = dACC_logFC,
        DLPFC_logFC = DLPFC_logFC
    )

    df_layer$color <- with(df_layer, ifelse(
        dACC_logFC > 1.25 & DLPFC_logFC > 1.25, "DLPFC & dACC enriched",
        ifelse(dACC_logFC > 1.25 & DLPFC_logFC < 1.25, "dACC enriched", "Neither")
    ))

    # genes that are either "DLPFC & dACC enriched" or "dACC enriched"
    write.csv(df_layer[df_layer$color != "Neither",],
              here("processed-data", "11_differential_expression", "novel_markers", paste0("dACC_DLPFC_", layer, ".csv")))

    if(layer == "L6a"){
        genes_to_label <- c("ISLR", "NR4A2", "DACH1", "KCTD8","TBR1")
    }
    if(layer == "L6b"){
        genes_to_label <- c("SEMA3A", "NXPH3", "ADRA2A", "SCUBE1", "CPLX3", "CRHBP")
    }

    # plot
    p <- ggplot(df_layer, aes(x = dACC_logFC, y = DLPFC_logFC)) +
        geom_point(aes(color = color), size=0.7) +
        scale_color_manual(values = c("DLPFC & dACC enriched" = "blue", "dACC enriched" = "red", "Neither" = "grey")) +
        geom_text_repel(aes(label = ifelse(gene %in% genes_to_label, gene, "")), size = 3, max.overlaps = Inf) +
        geom_smooth(method = "lm", se = FALSE, color = "blue") +
        geom_vline(xintercept = 1.25, color = "black", linetype = "dashed") +
        geom_hline(yintercept = 1.25, color = "black", linetype = "dashed") +
        geom_vline(xintercept = -1.25, color = "black", linetype = "dashed") +
        geom_hline(yintercept = -1.25, color = "black", linetype = "dashed") +
        xlim(-5, 5) +
        ylim(-5, 5) +
        labs(
            title = paste0("DLPFC: L6 & dACC Layer: ", layer),
            x = "dACC logFC",
            y = "DLPFC logFC",
            subtitle = paste0("Pearson's r: ", round(corr_layer$estimate, 2))
        )

    print(p)
    p_list[[layer]] <- p
}


dev.off()

pdf(file = here("plots", "11_differential_expression", "layer_markers_dACC_DLPFC.pdf"), width = 10, height = 10)
for (layer in names(p_list)) {
    print(p_list[[layer]])
}
dev.off()

