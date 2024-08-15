setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
# Load libraries
library(RcppML)
library(SingleCellExperiment)
library(here)
library(scuttle)
library(Matrix)
library(singlet)

# want to use logcounts for NMF
logcounts_combined <- readRDS(file = here("processed-data", "18_PsychENCODE_NMF", "pseudobulk", "pseudobulk_combined.rds"))

cvnmf <- cross_validate_nmf(
    logcounts_combined,
    ranks=c(25,50,75,100,125,150),
    n_replicates = 3,
    tol = 1e-03,
    maxit = 100,
    verbose = 3,
    L1 = 0.1,
    L2 = 0,
    threads = 0,
    test_density = 0.2
)

saveRDS(cvnmf, file = here("processed-data", "18_PsychENCODE_NMF", "nmf_cv_results.rds"))

pdf(here("plots", "18_PsychENCODE_NMF", "DLPFC_nmf_cv_results.pdf"))
plot(cvnmf)
dev.off()

