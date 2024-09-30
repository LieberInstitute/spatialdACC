setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("lobstr"))
suppressPackageStartupMessages(library("sessioninfo"))

load(here("processed-data", "VSPG","01_QC", "spe.RData"), verbose = TRUE)
dim(spe)
# 36601 19968
lobstr::obj_size(spe)
# 890.58 MB

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)
length(no_expr)
# [1] 8505
length(no_expr) / nrow(spe) * 100
# [1] 23.23707
spe <- spe[-no_expr, ]
dim(spe)
# [1] 28096 19968

## Now drop the spots outside the tissue
spe <- spe[, colData(spe)$in_tissue]
dim(spe)
# [1] 28096 17349

## Remove spots without counts
if (any(colSums(counts(spe)) == 0)) {
    message("removing spots without counts for spe")
    spe <- spe[, -which(colSums(counts(spe)) == 0)]
    dim(spe)
}

lobstr::obj_size(spe)
# 827.82 MB

save(spe, file = here::here("processed-data", "VSPG","01_QC", "spe_basic.Rdata"))
