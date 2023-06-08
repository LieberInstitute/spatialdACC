setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("ggspavis"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("here"))

load(here("processed-data", "05_QC", "spe_QC.Rdata"), verbose = TRUE)

## plot of tissue spots by percentage of reads mapped to mitochondrial transcripts
pdf(here("plots", "05_QC", "QC_mito_percent.pdf"), width = 21, height = 10)
samples <- unique(colData(spe)[, c("sample_id", "brnum")])
rownames(samples) <- NULL

for (i in 1:nrow(samples)) {
    p <- vis_gene(
        spe = spe,
        sampleid = samples$sample_id[i],
        geneid = "subsets_Mito_percent",
        spatial = T,
        point_size = 5,
        ... = paste0("_", samples$brnum[i])
    )

    p1 <- plotVisium(spe[, which(spe$sample_id == samples$sample_id[i])], spots = FALSE)

    grid.arrange(p, p1, nrow = 1)
}
dev.off()
