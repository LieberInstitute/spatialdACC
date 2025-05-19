setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("scater")
    library("spatialLIBD")
    library("dplyr")
    library("tidyr")
    library("ggplot2")
    library("patchwork")
})

nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", 9)
load(
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_pseudobulk", paste0(nnSVG_precast_name,".Rdata"))
)
mat_dACC <- assay(spe_pseudo, "logcounts")
groups_dACC <- factor(colData(spe_pseudo)[["layer"]])

sce <- fetch_data("sce_layer")
mat_dlPFC <- assay(sce, "logcounts")
groups_dlPFC <- factor(colData(sce)[["layer_guess_reordered_short"]])

#genes <- c("RELN", "MSX1", "VIM","HBB","NTS","HBA1")
#genes <- c("STXBP6", "LAMP5", "KCTD4", "ARHGAP4", "RSPO2", "CCNO","C1QL2")
#genes <- c("LINC01007", "ADCYAP1")
#genes <- c("RORB", "UNC5D", "PVALB")
#genes <- c("PCP4", "TRABD2A", "MEPE", "CD24", "CD52", "FDPS", "DRD5", "GYG1", "ITGB1BP1")
genes <- c("ISLR", "NR4A2", "DACH1", "KCTD8", "TBR1")
#genes <- c("SEMA3A", "NXPH3", "ADRA2A", "SCUBE1", "CPLX3", "CRHBP")

for(gene in genes){
    print(gene)
    pdf(file = here("plots", "11_differential_expression", "pseudobulk","boxplots_annotations", "L6a",
                    paste0(gene,".pdf")), height=10, width=15)

    par(mfrow = c(1, 2))  # Set plotting area to have 1 row and 2 columns

    i_dACC <- which(rowData(spe_pseudo)$gene_name == gene)
    i_dlPFC <- which(rowData(sce)$gene_name == gene)

    boxplot(
        mat_dACC[i_dACC, ] ~ groups_dACC,
        xlab = "",
        ylab = "",
        main = paste0(gene," expr. in dACC"),
        outline = FALSE,
        cex = 2,
        cex.axis = 2 * 4 / 5,
        cex.lab = 2,
        cex.main = ifelse(T, 2, 2 * 3 / 4),
        ylim = range(mat_dACC[i_dACC, ]),
        las = 2
    )
    points(
        mat_dACC[i_dACC, ] ~ jitter(as.integer(groups_dACC)),
        pch = 21,
        cex = 2
    )


    boxplot(
        mat_dlPFC[i_dlPFC, ] ~ groups_dlPFC,
        xlab = "",
        ylab = "",
        main = paste0(gene," expr. in dlPFC"),
        outline = FALSE,
        cex = 2,
        cex.axis = 2 * 4 / 5,
        cex.lab = 2,
        cex.main = ifelse(T, 2, 2 * 3 / 4),
        ylim = range(mat_dlPFC[i_dlPFC, ]),
        las = 2
    )
    points(
        mat_dlPFC[i_dlPFC, ] ~ jitter(as.integer(groups_dlPFC)),
        pch = 21,
        cex = 2
    )

    dev.off()
}





# subfigure
nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", 9)
load(
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_pseudobulk", paste0(nnSVG_precast_name,".Rdata"))
)
# we need to manually compute log2(cpm+1) instead of scaled pseudocount
assay(spe_pseudo, "logcounts") <- NULL
colData(spe_pseudo)$sizeFactor <- NULL
assay(spe_pseudo, "normounts") <- edgeR::cpm(counts(spe_pseudo))
assay(spe_pseudo, "logcounts") <- log2(assay(spe_pseudo, "normounts")+1)
mat_dACC <- assay(spe_pseudo, "logcounts")
groups_dACC <- factor(colData(spe_pseudo)[["layer"]])
spe_pseudo_dACC <- spe_pseudo

load(
    file = here("processed-data", "08_clustering", "DRD5_DLPFC_30_pseudobulk.Rdata")
)
assay(DRD5_DLPFC_30_pseudo, "normounts") <- edgeR::cpm(counts(DRD5_DLPFC_30_pseudo))
assay(DRD5_DLPFC_30_pseudo, "logcounts") <- log2(assay(DRD5_DLPFC_30_pseudo, "normounts")+1)
mat_dlPFC <- assay(DRD5_DLPFC_30_pseudo, "logcounts")
groups_dlPFC <- factor(colData(DRD5_DLPFC_30_pseudo)[["BayesSpace_harmony_09"]])
spe_pseudo_dlPFC <- DRD5_DLPFC_30_pseudo

genes <- c("CPLX3", "KCTD8", "DRD5")

dACC_layer_colors <- c(
    "L2" = "#377EB8",
    "L3" = "#4DAF4A",
    "L5" = "#FFD700",
    "L6b" = "#c46200",
    "L6a" = "#FFC18A",
    "WM" = "#1A1A1A",
    "L1" = "#F0027F"
)
dlPFC_layer_colors <- c(
    "L2" = "#377EB8",
    "L3" = "#4DAF4A",
    "L5" = "#FFD700",
    "L4" = "#984EA3",
    "L6" = "#FF7F00",
    "WM" = "#1A1A1A",
    "L1" = "#F0027F"
)

plot_list <- c()

for (gene in genes) {
    print(gene)
    i_dACC <- which(rowData(spe_pseudo_dACC)$gene_name == gene)
    i_dlPFC <- which(rowData(spe_pseudo_dlPFC)$gene_name == gene)

    df_dACC <- as.data.frame(mat_dACC[i_dACC, ])
    df_dACC$Layer <- groups_dACC
    colnames(df_dACC)[1] <- "Expression"

    df_dlPFC <- as.data.frame(mat_dlPFC[i_dlPFC, ])
    df_dlPFC$Layer <- groups_dlPFC
    colnames(df_dlPFC)[1] <- "Expression"

    p1 <- ggplot(df_dACC,aes(x=Layer, y=Expression, color=Layer)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(size = 0.6, alpha = 0.7) +
        theme_bw() +
        scale_color_manual(values=dACC_layer_colors) +
        ylab(expression("dACC " ~ log[2](cpm + 1))) +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        ggtitle(gene) +
        xlab("") +
        ylim(c(0,7.5))

    p2 <- ggplot(df_dlPFC,aes(x=Layer, y=Expression, color=Layer)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(size = 0.6, alpha = 0.7) +
        theme_bw() +
        scale_color_manual(values=dlPFC_layer_colors) +
        ylab(expression("dlPFC " ~ log[2](cpm + 1))) +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        xlab("") +
        ylim(c(0,7.5))


    plot_list <- c(plot_list, list(p1, p2))

}

pdf(file = here("plots", "11_differential_expression", "pseudobulk","boxplots_annotations",
                "dACC_dlPFC_pseudobulked.pdf"), height=5, width=7)

wrap_plots(plot_list[[1]], plot_list[[3]],
           plot_list[[5]], plot_list[[2]],
           plot_list[[4]], plot_list[[6]],
           nrow = 2, guides="collect") +
    plot_layout(axes = "collect_y")

dev.off()
