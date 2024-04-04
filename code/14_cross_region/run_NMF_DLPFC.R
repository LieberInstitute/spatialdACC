setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(RcppML)
library(pheatmap)
library(SingleCellExperiment)
library(here)
library(scuttle)
library(Matrix)

# Load data
spe <- spatialLIBD::fetch_data(type = "spe")

# want to use logcounts for NMF
spe <- logNormCounts(spe)

options(RcppML.threads=4)
x <- RcppML::nmf(assay(spe,'logcounts'),
                 tol = 1e-06,
                 maxit = 1000,
                 verbose = T,
                 L1 = 0.1,
                 seed = 1135,
                 mask_zeros = FALSE,
                 diag = TRUE,
                 nonneg = TRUE
)

# Save NMF results
saveRDS(x, file = here("processed-data", "14_cross_region", "DLPFC_nmf_results.RDS"))


