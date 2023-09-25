setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("BayesSpace")
    library("RColorBrewer")
    library("ggplot2")
    library("gridExtra")
    library("Polychrome")
})

load(here("processed-data", "07_batch_correction", "spe_harmony.Rdata"))
spe_harmony <- spe
dim(spe_harmony)

spe_harmony <- qTune(spe_harmony, qs=seq(5, 20), platform="ST")

load(here("processed-data", "07_batch_correction", "spe_mnn.Rdata"))
spe_mnn <- spe
dim(spe_mnn)

spe_mnn <- qTune(spe_mnn, qs=seq(5, 20), platform="ST")

pdf(file = here::here("plots", "08_clustering", "BayesSpace", "qtune.pdf"), width = 21, height = 20)
qPlot(spe_harmony)
qPlot(spe_mnn)
dev.off()
