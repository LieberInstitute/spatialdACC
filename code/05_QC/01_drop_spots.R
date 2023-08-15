setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("lobstr"))
suppressPackageStartupMessages(library("sessioninfo"))

load(here("processed-data", "02_build_spe", "spe_raw.Rdata"), verbose = TRUE)
dim(spe)
# 36601 84864
lobstr::obj_size(spe)
# 3.33 GB

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)
length(no_expr)
# [1] 6881
length(no_expr) / nrow(spe) * 100
# [1] 18.80003
spe <- spe[-no_expr, ]
dim(spe)
# [1] 29720 84864

## Now drop the spots outside the tissue
spe <- spe[, colData(spe)$in_tissue]
dim(spe)
# [1] 29720 77604

## Remove spots without counts
if (any(colSums(counts(spe)) == 0)) {
    message("removing spots without counts for spe")
    spe <- spe[, -which(colSums(counts(spe)) == 0)]
    dim(spe)
}

# [1] 29720 77604

lobstr::obj_size(spe)
# 3.29 GB

save(spe, file = here::here("processed-data", "05_QC", "spe_basic.Rdata"))
