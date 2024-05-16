setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(scater)
library(scran)

# based on
# https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/misc.html

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))

sce <- computeSumFactors(sce)
sce <- logNormCounts(sce)

correlatePairs(sce, subset.row=cbind("VAT1L", "DRD5"))

pdf(file = here::here("plots", "16_VENs_analysis", "single_nucleus_coexpression.pdf"),
    width = 21, height = 20)
plotExpression(sce, features="VAT1L", x="DRD5")
dev.off()
