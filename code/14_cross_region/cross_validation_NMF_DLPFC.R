setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(RcppML)
library(SingleCellExperiment)
library(here)
library(scuttle)
library(Matrix)
library(singlet)
library(spatialLIBD)

# Load data
spe <- spatialLIBD::fetch_data(type = "spe")

# want to use logcounts for NMF
spe <- logNormCounts(spe)

cvnmf <- cross_validate_nmf(
    logcounts(spe),
    ranks=c(15,25,35,45,55,75,100),
    n_replicates = 3,
    tol = 1e-03,
    maxit = 100,
    verbose = 3,
    L1 = 0.1,
    L2 = 0,
    threads = 0,
    test_density = 0.2
)

saveRDS(cvnmf, file = here("processed-data", "14_cross_region", "DLPFC_nmf_cv_results.RDS"))

pdf(here("plots", "14_cross_region", "DLPFC_nmf_cv_results.pdf"))
plot(cvnmf)
dev.off()

