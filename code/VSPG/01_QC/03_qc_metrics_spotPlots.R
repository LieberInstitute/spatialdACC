setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("ggspavis"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("here"))

load(here("processed-data", "VSPG", "01_QC", "spe_QC.Rdata"), verbose = TRUE)

## QC plot of tissue spots discarded by brain with low sum and features
pdf(here("plots", "VSPG", "01_QC", "QC_brain.pdf"), width = 21, height = 20)
brains <- unique(spe$brnum)

for (i in seq_along(brains)){
    speb <- spe[, which(spe$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    print(length(samples))

    p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 2,... = paste0("_",brains[i]) ) +
        ggtitle("all metrics combined")
    p2 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "high_mito_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 2,... = paste0("_",brains[i]) ) +
        ggtitle("high mito")
    p3 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "low_sum_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 2,... = paste0("_",brains[i]) ) +
        ggtitle("low sum")
    p4 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "low_detected_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 2,... = paste0("_",brains[i]) ) +
        ggtitle("low detected")
    p <- plotVisium(spe[, which(spe$sample_id == samples)], spots = FALSE)

    grid.arrange(p, p1, p2, p3, p4, nrow = 3)

}

dev.off()
