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
    library("patchwork")
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

# old plots code
pdf(file = here("plots", "11_differential_expression", "layer_markers_dACC_DLPFC.pdf"), width = 7, height = 7)

# repeat this in a loop for layers L1, L2, L3, L5
layers <- c("L1","L2", "L3", "L5")
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
        dACC_logFC > 1.25 & DLPFC_logFC > 1.25, "dlPFC & dACC enriched",
        ifelse(dACC_logFC > 1.25 & DLPFC_logFC < 1.25, "dACC enriched",
               ifelse(DLPFC_logFC > 1.25 & dACC_logFC < 1.25, "dlPFC enriched", "None")
        )
    ))


    # write to csv file
    # genes that are either "DLPFC & dACC enriched" or "dACC enriched"
    write.csv(df_layer[df_layer$color != "None",],
              here("processed-data", "11_differential_expression", "novel_markers", paste0("dACC_DLPFC_", layer, ".csv")))

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
        genes_to_label <- c("PCP4", "TRABD2A", "MEPE", "CD24", "CD52", "FDPS", "DRD5", "GYG1", "ITGB1BP1",
                            "VAT1L", "SULF2", "HAPLN4", "GABRQ", "FEZF2", "POU3F1")
    }

    # plot
    p <- ggplot(df_layer, aes(x = dACC_logFC, y = DLPFC_logFC)) +
        geom_point(aes(color = color), size=0.5) +
        scale_color_manual(values = c("dlPFC & dACC enriched" = "blue", "dACC enriched" = "red", "Neither" = "grey")) +
        geom_text_repel(aes(label = ifelse(gene %in% genes_to_label, gene, "")), size = 3, max.overlaps = Inf) +
        geom_smooth(method = "lm", se = FALSE, color = "blue") +
        geom_vline(xintercept = 1.25, color = "black", linetype = "dashed") +
        geom_hline(yintercept = 1.25, color = "black", linetype = "dashed") +
        geom_vline(xintercept = -1.25, color = "black", linetype = "dashed") +
        geom_hline(yintercept = -1.25, color = "black", linetype = "dashed") +
        xlim(-3.5, 3.5) +
        ylim(-3.5, 3.5) +
        labs(
            title = paste0("Layer: ", layer),
            x = "dACC logFC",
            y = "dlPFC logFC",
            subtitle = paste0("Pearson's r: ", round(corr_layer$estimate, 2))
        ) +
        theme_bw()

    #print(p)
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
        dACC_logFC > 1.25 & DLPFC_logFC > 1.25, "dlPFC & dACC enriched",
        ifelse(dACC_logFC > 1.25 & DLPFC_logFC < 1.25, "dACC enriched",
               ifelse(DLPFC_logFC > 1.25 & dACC_logFC < 1.25, "dlPFC enriched", "None")
        )
    ))


    # genes that are either "DLPFC & dACC enriched" or "dACC enriched"
    write.csv(df_layer[df_layer$color != "None",],
              here("processed-data", "11_differential_expression", "novel_markers", paste0("dACC_DLPFC_", layer, ".csv")))

    if(layer == "L6a"){
        genes_to_label <- c("ISLR", "NR4A2", "DACH1", "KCTD8","TBR1")
    }
    if(layer == "L6b"){
        genes_to_label <- c("SEMA3A", "NXPH3", "ADRA2A", "SCUBE1", "CPLX3", "CRHBP")
    }

    # plot
    p <- ggplot(df_layer, aes(x = dACC_logFC, y = DLPFC_logFC)) +
        geom_point(aes(color = color), size=0.5) +
        scale_color_manual(values = c("dlPFC & dACC enriched" = "blue", "dACC enriched" = "red", "Neither" = "grey")) +
        geom_text_repel(aes(label = ifelse(gene %in% genes_to_label, gene, "")), size = 3, max.overlaps = Inf) +
        geom_smooth(method = "lm", se = FALSE, color = "blue") +
        geom_vline(xintercept = 1.25, color = "black", linetype = "dashed") +
        geom_hline(yintercept = 1.25, color = "black", linetype = "dashed") +
        geom_vline(xintercept = -1.25, color = "black", linetype = "dashed") +
        geom_hline(yintercept = -1.25, color = "black", linetype = "dashed") +
        xlim(-3.5, 3.5) +
        ylim(-3.5, 3.5) +
        labs(
            title = paste0("dlPFC: L6 & dACC Layer: ", layer),
            x = "dACC logFC",
            y = "dlPFC logFC",
            subtitle = paste0("Pearson's r: ", round(corr_layer$estimate, 2))
        ) +
        theme_bw()

    #print(p)
    p_list[[layer]] <- p
}


dev.off()











# final plots code
layer <- "L5"
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
    dACC_logFC > 1.25 & DLPFC_logFC > 1.25, "dlPFC & dACC enriched",
    ifelse(dACC_logFC > 1.25 & DLPFC_logFC < 1.25, "dACC enriched",
           ifelse(DLPFC_logFC > 1.25 & dACC_logFC < 1.25, "dlPFC enriched", "None")
    )
))

genes_to_label_L5 <- c("PCP4", "TRABD2A", "MEPE", "CD24", "CD52", "FDPS", "DRD5", "GYG1", "ITGB1BP1",
                        "VAT1L", "SULF2", "HAPLN4", "GABRQ", "FEZF2", "POU3F1")

# plot
p5 <- ggplot(df_layer, aes(x = dACC_logFC, y = DLPFC_logFC)) +
    geom_point(aes(color = color), size=0.5) +
    scale_color_manual(values = c("dlPFC & dACC enriched" = "purple",
                                  "dACC enriched" = "red",
                                  "dlPFC enriched" = "blue",
                                  "Neither" = "grey")) +
    geom_text_repel(aes(label = ifelse(gene %in% genes_to_label_L5, gene, "")), size = 3, max.overlaps = Inf) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    geom_vline(xintercept = 1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 1.25, color = "black", linetype = "dashed") +
    geom_vline(xintercept = -1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = -1.25, color = "black", linetype = "dashed") +
    xlim(-2.6, 3.9) +
    ylim(-1.7, 2.7) +
    labs(
        title = "dlPFC L5 & dACC L5",
        x = "dACC logFC",
        y = "dlPFC logFC",
        subtitle = paste0("Pearson's r: ", round(corr_layer$estimate, 2))
    ) +
    theme_bw()

layer <- "L6a"
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
    dACC_logFC > 1.25 & DLPFC_logFC > 1.25, "dlPFC & dACC enriched",
    ifelse(dACC_logFC > 1.25 & DLPFC_logFC < 1.25, "dACC enriched",
           ifelse(DLPFC_logFC > 1.25 & dACC_logFC < 1.25, "dlPFC enriched", "None")
    )
))
genes_to_label_L6a <- c("ISLR", "NR4A2", "DACH1", "KCTD8","TBR1")

# plot
p6a <- ggplot(df_layer, aes(x = dACC_logFC, y = DLPFC_logFC)) +
    geom_point(aes(color = color), size=0.5) +
    scale_color_manual(values = c("dlPFC & dACC enriched" = "purple",
                                  "dACC enriched" = "red",
                                  "dlPFC enriched" = "blue",
                                  "Neither" = "grey")) +
    geom_text_repel(aes(label = ifelse(gene %in% genes_to_label_L6a, gene, "")), size = 3, max.overlaps = Inf) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    geom_vline(xintercept = 1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 1.25, color = "black", linetype = "dashed") +
    geom_vline(xintercept = -1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = -1.25, color = "black", linetype = "dashed") +
    xlim(-2.6, 3.9) +
    ylim(-1.7, 2.7) +
    labs(
        title = paste0("dlPFC L6 & dACC ", layer),
        x = "dACC logFC",
        y = "dlPFC logFC",
        subtitle = paste0("Pearson's r: ", round(corr_layer$estimate, 2))
    ) +
    theme_bw()

layer <- "L6b"
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
    dACC_logFC > 1.25 & DLPFC_logFC > 1.25, "dlPFC & dACC enriched",
    ifelse(dACC_logFC > 1.25 & DLPFC_logFC < 1.25, "dACC enriched",
           ifelse(DLPFC_logFC > 1.25 & dACC_logFC < 1.25, "dlPFC enriched", "None")
    )
))
genes_to_label_L6b <- c("SEMA3A", "NXPH3", "ADRA2A", "SCUBE1", "CPLX3", "CRHBP")

# plot
p6b <- ggplot(df_layer, aes(x = dACC_logFC, y = DLPFC_logFC)) +
    geom_point(aes(color = color), size=0.5) +
    scale_color_manual(values = c("dlPFC & dACC enriched" = "purple",
                                  "dACC enriched" = "red",
                                  "dlPFC enriched" = "blue",
                                  "Neither" = "grey")) +
    geom_text_repel(aes(label = ifelse(gene %in% genes_to_label_L6b, gene, "")), size = 3, max.overlaps = Inf) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    geom_vline(xintercept = 1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 1.25, color = "black", linetype = "dashed") +
    geom_vline(xintercept = -1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = -1.25, color = "black", linetype = "dashed") +
    xlim(-2.6, 3.9) +
    ylim(-1.7, 2.7) +
    labs(
        title = paste0("dlPFC L6 & dACC ", layer),
        x = "dACC logFC",
        y = "dlPFC logFC",
        subtitle = paste0("Pearson's r: ", round(corr_layer$estimate, 2))
    ) +
    theme_bw()

pdf(file = here("plots", "11_differential_expression", "layer_markers_dACC_DLPFC_subset.pdf"), width = 14, height = 7)
wrap_plots(p5, p6a, p6b, nrow=1, guides = "collect") + plot_layout(axes = "collect") & theme(legend.position = "bottom")
dev.off()

