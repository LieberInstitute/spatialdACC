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

genes <- c("VAT1L")

for(gene in genes){
    print(gene)
    png(file = here("plots", "11_differential_expression", "pseudobulk","boxplots_annotations",
                    paste0(gene,".png")), height=8, width=8, unit="in",res=300)


    i_dACC <- which(rowData(spe_pseudo)$gene_name == gene)

    boxplot(
        mat_dACC[i_dACC, ] ~ groups_dACC,
        xlab = "",
        ylab = "",
        main = "",
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

    dev.off()
}
