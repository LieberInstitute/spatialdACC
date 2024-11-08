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
#0.12

pdf(file = here::here("plots", "16_VENs_analysis", "single_nucleus_coexpression.pdf"),
    width = 21, height = 20)
plotExpression(sce, features="VAT1L", x="DRD5")
dev.off()


correlatePairs(sce, subset.row=cbind("VAT1L", "GABRQ"))
#0.14
pdf(file = here::here("plots", "16_VENs_analysis", "single_nucleus_coexpression_VAT1L_GABRQ.pdf"),
    width = 21, height = 20)
plotExpression(sce, features="VAT1L", x="GABRQ")
dev.off()

correlatePairs(sce, subset.row=cbind("VAT1L", "POU3F1"))
#0.11
correlatePairs(sce, subset.row=cbind("VAT1L", "ADRA1A"))
#0.23
