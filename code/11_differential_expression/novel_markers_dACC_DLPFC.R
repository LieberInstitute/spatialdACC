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

load(file = here("processed-data", "11_differential_expression",paste0("DLPFC_30_DE_with_DRD5",".Rdata")))
modeling_results_DLPFC <- modeling_results

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

genes_to_label_L5 <- c("RORB", "PVALB", "HAPLN4", "SULF2", "POU3F1", "GRIN3A", "OPRM1",
                       "HTR2C", "FOXP2", "TSHZ2", "DRD5", "PCP4", "FEZF2", "TRABD2A", "MEPE")

# plot
p5 <- ggplot(df_layer, aes(x = dACC_logFC, y = DLPFC_logFC)) +
    geom_point(aes(color = color), size=0.5) +
    scale_color_manual(values = c("dlPFC & dACC enriched" = "purple",
                                  "dACC enriched" = "red",
                                  "dlPFC enriched" = "blue",
                                  "Neither" = "grey")) +
    geom_label_repel(
        aes(label = ifelse(gene %in% genes_to_label_L5, gene, "")),
        size = 2,
        max.overlaps = Inf,
        force = 2,
        force_pull = 0.5,
        box.padding = 0.5,
        point.padding = 0.3,
        max.iter = 10000,
        nudge_y = 0.1,
        nudge_x = -0.1
    ) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    geom_vline(xintercept = 1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 1.25, color = "black", linetype = "dashed") +
    geom_vline(xintercept = -1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = -1.25, color = "black", linetype = "dashed") +
    xlim(-2.6, 3.9) +
    ylim(-1.7, 2.9) +
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
genes_to_label_L6a <- c("KCTD8", "TBR1", "SSTR2", "SEMA3A", "NR4A2", "DACH1", "ISLR",
                        "KRT17", "SEMA3E", "NXPH4")

# plot
p6a <- ggplot(df_layer, aes(x = dACC_logFC, y = DLPFC_logFC)) +
    geom_point(aes(color = color), size=0.5) +
    scale_color_manual(values = c("dlPFC & dACC enriched" = "purple",
                                  "dACC enriched" = "red",
                                  "dlPFC enriched" = "blue",
                                  "Neither" = "grey")) +
    geom_label_repel(
        aes(label = ifelse(gene %in% genes_to_label_L6a, gene, "")),
        size = 2,
        max.overlaps = Inf,
        force = 2,
        force_pull = 0.5,
        box.padding = 0.5,
        point.padding = 0.3,
        max.iter = 10000,
        nudge_y = 0.1
    ) +
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
genes_to_label_L6b <- c("CPLX3", "CRHBP", "ADRA2A", "NOS1", "ISLR", "KRT17", "SEMA3A",
                        "FOXP2", "NXPH4")

# plot
p6b <- ggplot(df_layer, aes(x = dACC_logFC, y = DLPFC_logFC)) +
    geom_point(aes(color = color), size=0.5) +
    scale_color_manual(values = c("dlPFC & dACC enriched" = "purple",
                                  "dACC enriched" = "red",
                                  "dlPFC enriched" = "blue",
                                  "Neither" = "grey")) +
    geom_label_repel(
        aes(label = ifelse(gene %in% genes_to_label_L6b, gene, "")),
        size = 2,
        max.overlaps = Inf,
        force = 2,
        force_pull = 0.5,
        box.padding = 0.5,
        point.padding = 0.3,
        max.iter = 10000
    ) +
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

layer <- "L2"
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

genes_to_label_L2 <- c("GULP1", "RSPO2", "BDNF", "SYTL5", "SHISA8", "CARTPT",
                       "HPCAL1", "LAMP5", "C1QL2", "CUX2", "FREM3")
# plot
p2 <- ggplot(df_layer, aes(x = dACC_logFC, y = DLPFC_logFC)) +
    geom_point(aes(color = color), size=0.5) +
    scale_color_manual(values = c("dlPFC & dACC enriched" = "purple",
                                  "dACC enriched" = "red",
                                  "dlPFC enriched" = "blue",
                                  "Neither" = "grey")) +
    geom_label_repel(
        aes(label = ifelse(gene %in% genes_to_label_L2, gene, "")),
        size = 2,
        max.overlaps = Inf,
        force = 2,
        force_pull = 0.5,
        box.padding = 0.5,
        point.padding = 0.3,
        max.iter = 10000,
        nudge_y = 0.2,
        nudge_x = -0.05
    ) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    geom_vline(xintercept = 1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 1.25, color = "black", linetype = "dashed") +
    geom_vline(xintercept = -1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = -1.25, color = "black", linetype = "dashed") +
    labs(
        title = "dlPFC L2 & dACC L2",
        x = "dACC logFC",
        y = "dlPFC logFC",
        subtitle = paste0("Pearson's r: ", round(corr_layer$estimate, 2))
    ) +
    theme_bw()

pdf(file = here("plots", "11_differential_expression", "layer_markers_dACC_DLPFC_subset.pdf"), width = 13, height = 5)
wrap_plots(p2, p5, p6a, p6b, nrow=1, guides = "collect") + plot_layout(axes = "collect") & theme(legend.position = "bottom")
dev.off()


# supp figure
layer <- "L1"
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

genes_to_label_L1 <- c("RELN", "MSX1", "VIM","HBB","NTS","HBA1")

# plot
p1 <- ggplot(df_layer, aes(x = dACC_logFC, y = DLPFC_logFC)) +
    geom_point(aes(color = color), size=0.5) +
    scale_color_manual(values = c("dlPFC & dACC enriched" = "purple",
                                  "dACC enriched" = "red",
                                  "dlPFC enriched" = "blue",
                                  "Neither" = "grey")) +
    geom_text_repel(aes(label = ifelse(gene %in% genes_to_label_L1, gene, "")), size = 3, max.overlaps = Inf) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    geom_vline(xintercept = 1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 1.25, color = "black", linetype = "dashed") +
    geom_vline(xintercept = -1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = -1.25, color = "black", linetype = "dashed") +
    labs(
        title = "dlPFC L1 & dACC L1",
        x = "dACC logFC",
        y = "dlPFC logFC",
        subtitle = paste0("Pearson's r: ", round(corr_layer$estimate, 2))
    ) +
    theme_bw()

layer <- "L3"
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

genes_to_label_L3 <- c("LINC01007", "ADCYAP1")

p3 <- ggplot(df_layer, aes(x = dACC_logFC, y = DLPFC_logFC)) +
    geom_point(aes(color = color), size=0.5) +
    scale_color_manual(values = c("dlPFC & dACC enriched" = "purple",
                                  "dACC enriched" = "red",
                                  "dlPFC enriched" = "blue",
                                  "Neither" = "grey")) +
    geom_text_repel(aes(label = ifelse(gene %in% genes_to_label_L3, gene, "")), size = 3, max.overlaps = Inf) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    geom_vline(xintercept = 1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 1.25, color = "black", linetype = "dashed") +
    geom_vline(xintercept = -1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = -1.25, color = "black", linetype = "dashed") +
    labs(
        title = "dlPFC L3 & dACC L3",
        x = "dACC logFC",
        y = "dlPFC logFC",
        subtitle = paste0("Pearson's r: ", round(corr_layer$estimate, 2))
    ) +
    theme_bw() +
    theme(legend.position="none")

png(file = here("plots", "11_differential_expression", "layer_markers_dACC_DLPFC_other.png"), width = 9, height = 5, unit="in",res=300)
wrap_plots(p1, p3, nrow=1, guides = "collect") + plot_layout(axes = "collect") & theme(legend.position = "bottom")
dev.off()





# supp figure 2
layer_dACC <- "L3"
layer_dlPFC <- "L4"
genes <- dACC_results[, "gene"]
dACC_logFC <- dACC_results[, paste0("logFC_", layer_dACC)]
DLPFC_logFC <- DLPFC_results[, paste0("logFC_", layer_dlPFC)]
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

genes_to_label <- c("RORB", "PVALB", "HAPLN4", "SULF2", "POU3F1", "GRIN3A", "OPRM1",
                    "HTR2C", "FOXP2", "TSHZ2", "DRD5", "PCP4", "FEZF2", "TRABD2A",
                    "MEPE", "LINC01007", "ADCYAP1")

# plot
p1 <- ggplot(df_layer, aes(x = dACC_logFC, y = DLPFC_logFC)) +
    geom_point(aes(color = color), size=0.5) +
    scale_color_manual(values = c("dlPFC & dACC enriched" = "purple",
                                  "dACC enriched" = "red",
                                  "dlPFC enriched" = "blue",
                                  "Neither" = "grey")) +
    geom_text_repel(aes(label = ifelse(gene %in% genes_to_label, gene, "")), size = 3, max.overlaps = Inf) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    geom_vline(xintercept = 1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 1.25, color = "black", linetype = "dashed") +
    geom_vline(xintercept = -1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = -1.25, color = "black", linetype = "dashed") +
    labs(
        title = "dlPFC L4 & dACC L3",
        x = "dACC logFC",
        y = "dlPFC logFC",
        subtitle = paste0("Pearson's r: ", round(corr_layer$estimate, 2))
    ) +
    theme_bw()

layer_dACC <- "L5"
layer_dlPFC <- "L4"
genes <- dACC_results[, "gene"]
dACC_logFC <- dACC_results[, paste0("logFC_", layer_dACC)]
DLPFC_logFC <- DLPFC_results[, paste0("logFC_", layer_dlPFC)]
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

genes_to_label <- c("RORB", "PVALB", "HAPLN4", "SULF2", "POU3F1", "GRIN3A", "OPRM1",
                    "HTR2C", "FOXP2", "TSHZ2", "DRD5", "PCP4", "FEZF2", "TRABD2A",
                    "MEPE", "LINC01007", "ADCYAP1")

p3 <- ggplot(df_layer, aes(x = dACC_logFC, y = DLPFC_logFC)) +
    geom_point(aes(color = color), size=0.5) +
    scale_color_manual(values = c("dlPFC & dACC enriched" = "purple",
                                  "dACC enriched" = "red",
                                  "dlPFC enriched" = "blue",
                                  "Neither" = "grey")) +
    geom_text_repel(aes(label = ifelse(gene %in% genes_to_label, gene, "")), size = 3, max.overlaps = Inf) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    geom_vline(xintercept = 1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 1.25, color = "black", linetype = "dashed") +
    geom_vline(xintercept = -1.25, color = "black", linetype = "dashed") +
    geom_hline(yintercept = -1.25, color = "black", linetype = "dashed") +
    labs(
        title = "dlPFC L4 & dACC L5",
        x = "dACC logFC",
        y = "dlPFC logFC",
        subtitle = paste0("Pearson's r: ", round(corr_layer$estimate, 2))
    ) +
    theme_bw() +
    theme(legend.position="none")

png(file = here("plots", "11_differential_expression", "layer_markers_dACC_DLPFC_L4.png"), width = 9, height = 5, unit="in",res=300)
wrap_plots(p1, p3, nrow=1, guides = "collect") + plot_layout(axes = "collect") & theme(legend.position = "bottom")
dev.off()
