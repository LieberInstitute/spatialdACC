setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("scater")
    library("spatialLIBD")
    library("dplyr")
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
genes <- c("PCP4", "TRABD2A", "MEPE", "CD24", "CD52", "FDPS", "DRD5", "GYG1", "ITGB1BP1")
#genes <- c("ISLR", "NR4A2", "DACH1", "KCTD8")
#genes <- c("SEMA3A", "NXPH3", "ADRA2A", "SCUBE1", "CPLX3", "CRHBP")

for(gene in genes){
    print(gene)
    pdf(file = here("plots", "11_differential_expression", "pseudobulk","boxplots_annotations", "L5",
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
