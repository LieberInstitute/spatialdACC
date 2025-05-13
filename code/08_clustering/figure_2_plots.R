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
    library("escheR")
})

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

load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata"))
spe_dACC <- spe

spe_DLPFC_30 <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium")
# create spatial labels for DLPFC_30
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 3] <- "L2"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 8] <- "L4"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 7] <- "L6"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 5] <- "L3"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 6] <- "WM"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 4] <- "L5"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 2] <- "L1"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 1] <- "meninges"
spe_DLPFC_30$BayesSpace_harmony_09[spe_DLPFC_30$BayesSpace_harmony_09 == 9] <- "WM"

# remove meninges
spe_DLPFC_30 <- spe_DLPFC_30[,which(spe_DLPFC_30$BayesSpace_harmony_09 != "meninges")]

spe_DLPFC_30$cluster <- spe_DLPFC_30$BayesSpace_harmony_09


# need to manually re-calculate to get DRD5 pseudobulked samples which were dropped due to filtering before
DRD5_DLPFC_30_pseudo <- scuttle::aggregateAcrossCells(
    spe_DLPFC_30,
    DataFrame(
        registration_variable = spe_DLPFC_30[["BayesSpace_harmony_09"]],
        registration_sample_id = spe_DLPFC_30[["sample_id"]]
    )
)
colnames(DRD5_DLPFC_30_pseudo) <-
    paste0(
        DRD5_DLPFC_30_pseudo$"sample_id",
        "_",
        DRD5_DLPFC_30_pseudo$"BayesSpace_harmony_09"
    )

colData(DRD5_DLPFC_30_pseudo)$sizeFactor <- NULL
assay(DRD5_DLPFC_30_pseudo, "normounts") <- edgeR::cpm(counts(DRD5_DLPFC_30_pseudo))
assay(DRD5_DLPFC_30_pseudo, "logcounts") <- log2(assay(DRD5_DLPFC_30_pseudo, "normounts")+1)
mat_dlPFC <- assay(DRD5_DLPFC_30_pseudo, "logcounts")
groups_dlPFC <- factor(colData(DRD5_DLPFC_30_pseudo)[["BayesSpace_harmony_09"]])


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

genes_exp_1 <- c("PCP4", "ADCYAP1", "RORB", "ARHGAP4")

boxplot_list_exp1 <- c()

for (gene in genes_exp_1) {
    print(gene)
    i_dACC <- which(rowData(spe_pseudo_dACC)$gene_name == gene)
    i_dlPFC <- which(rowData(DRD5_DLPFC_30_pseudo)$gene_name == gene)

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
        ylim(c(0,max(df_dACC$Expression))) +
        theme(legend.position="none")

    p2 <- ggplot(df_dlPFC,aes(x=Layer, y=Expression, color=Layer)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(size = 0.6, alpha = 0.7) +
        theme_bw() +
        scale_color_manual(values=dlPFC_layer_colors) +
        ylab(expression("dlPFC " ~ log[2](cpm + 1))) +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        ggtitle(gene) +
        xlab("") +
        ylim(c(0,max(df_dlPFC$Expression))) +
        theme(legend.position="none")


    boxplot_list_exp1 <- c(boxplot_list_exp1, list(p1, p2))

}

spotplot_list_exp1 <- c()

for (gene in genes_exp_1) {
    print(gene)

    spe_dACC_1 <- spe_dACC[, which(spe_dACC$sample_id == "V12N28-331_B1")]
    spe_dACC_1$counts_gene <- logcounts(spe_dACC_1)[which(rowData(spe_dACC_1)$gene_name==gene),]
    p1 <- make_escheR(spe_dACC_1) |> add_fill(var="counts_gene", point_size = 1) |> add_ground(var="layer", stroke=0.1, point_size = 1) +
        scale_fill_gradient(low = "white", high = "black") + labs(title = "") +
        scale_color_manual(values=dACC_layer_colors) +
        theme(legend.position="none")

    spe_DLPFC_30_1 <- spe_DLPFC_30[, which(spe_DLPFC_30$sample_id == "Br6432_ant")]
    spe_DLPFC_30_1$counts_gene <- logcounts(spe_DLPFC_30_1)[which(rowData(spe_DLPFC_30_1)$gene_name==gene),]
    p2 <- make_escheR(spe_DLPFC_30_1) |> add_fill(var="counts_gene", point_size = 1) |> add_ground(var="cluster", stroke=0.1, point_size = 1) +
        scale_fill_gradient(low = "white", high = "black") + labs(title = "") +
        scale_color_manual(values=dlPFC_layer_colors) +
        theme(legend.position="none")

    spotplot_list_exp1 <- c(spotplot_list_exp1, list(p1, p2))

}


genes_exp_2 <- c("CPLX3", "KCTD8", "NXPH3", "DRD5")

boxplot_list_exp2 <- c()

for (gene in genes_exp_2) {
    print(gene)
    i_dACC <- which(rowData(spe_pseudo_dACC)$gene_name == gene)
    i_dlPFC <- which(rowData(DRD5_DLPFC_30_pseudo)$gene_name == gene)

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
        ylim(c(0,max(df_dACC$Expression))) +
        theme(legend.position="none")

    p2 <- ggplot(df_dlPFC,aes(x=Layer, y=Expression, color=Layer)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(size = 0.6, alpha = 0.7) +
        theme_bw() +
        scale_color_manual(values=dlPFC_layer_colors) +
        ylab(expression("dlPFC " ~ log[2](cpm + 1))) +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        ggtitle(gene) +
        xlab("") +
        ylim(c(0,max(df_dlPFC$Expression))) +
        theme(legend.position="none")


    boxplot_list_exp2 <- c(boxplot_list_exp2, list(p1, p2))

}

spotplot_list_exp2 <- c()

for (gene in genes_exp_2) {
    print(gene)

    spe_dACC_1 <- spe_dACC[, which(spe_dACC$sample_id == "V12N28-331_B1")]
    spe_dACC_1$counts_gene <- logcounts(spe_dACC_1)[which(rowData(spe_dACC_1)$gene_name==gene),]
    p1 <- make_escheR(spe_dACC_1) |> add_fill(var="counts_gene", point_size = 1) |> add_ground(var="layer", stroke=0.1, point_size = 1) +
        scale_fill_gradient(low = "white", high = "black") + labs(title = "") +
        scale_color_manual(values=dACC_layer_colors) +
        theme(legend.position="none")

    spe_DLPFC_30_1 <- spe_DLPFC_30[, which(spe_DLPFC_30$sample_id == "Br6432_ant")]
    spe_DLPFC_30_1$counts_gene <- logcounts(spe_DLPFC_30_1)[which(rowData(spe_DLPFC_30_1)$gene_name==gene),]
    p2 <- make_escheR(spe_DLPFC_30_1) |> add_fill(var="counts_gene", point_size = 1) |> add_ground(var="cluster", stroke=0.1, point_size = 1) +
        scale_fill_gradient(low = "white", high = "black") + labs(title = "") +
        scale_color_manual(values=dlPFC_layer_colors) +
        theme(legend.position="none")

    spotplot_list_exp2 <- c(spotplot_list_exp2, list(p1, p2))

}




pdf(file = here("plots", "08_clustering",
                "dACC_dlPFC_figure_2.pdf"), height=13, width=6)

wrap_plots(boxplot_list_exp1[[1]], spotplot_list_exp1[[1]],
           boxplot_list_exp1[[3]], spotplot_list_exp1[[3]],
           boxplot_list_exp1[[5]], spotplot_list_exp1[[5]],
           boxplot_list_exp1[[7]], spotplot_list_exp1[[7]],
           nrow = 4)

wrap_plots(boxplot_list_exp1[[2]], spotplot_list_exp1[[2]],
           boxplot_list_exp1[[4]], spotplot_list_exp1[[4]],
           boxplot_list_exp1[[6]], spotplot_list_exp1[[6]],
           boxplot_list_exp1[[8]], spotplot_list_exp1[[8]],
           nrow = 4)

wrap_plots(boxplot_list_exp2[[1]], spotplot_list_exp2[[1]],
           boxplot_list_exp2[[3]], spotplot_list_exp2[[3]],
           boxplot_list_exp2[[5]], spotplot_list_exp2[[5]],
           boxplot_list_exp2[[7]], spotplot_list_exp2[[7]],
           nrow = 4)

wrap_plots(boxplot_list_exp2[[2]], spotplot_list_exp2[[2]],
           boxplot_list_exp2[[4]], spotplot_list_exp2[[4]],
           boxplot_list_exp2[[6]], spotplot_list_exp2[[6]],
           boxplot_list_exp2[[8]], spotplot_list_exp2[[8]],
           nrow = 4)


dev.off()
