setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
# Load libraries
library(RcppML)
library(SingleCellExperiment)
library(here)
library(scuttle)
library(Matrix)
library(singlet)

# want to use logcounts for NMF
sce <- readRDS(file = here("processed-data", "18_PsychENCODE_NMF", "pseudobulk_combined.rds"))

set.seed(12)
cvnmf <- cross_validate_nmf(
    logcounts(sce),
    ranks=c(25,50,75,100),
    n_replicates = 1,
    tol = 1e-03,
    maxit = 100,
    verbose = 3,
    L1 = 0.1,
    L2 = 0,
    threads = 0,
    test_density = 0.2
)

pdf(here("plots", "18_PsychENCODE_NMF", "nmf_cv_results.pdf"))
plot(cvnmf)
dev.off()

saveRDS(cvnmf, file = here("processed-data", "18_PsychENCODE_NMF", "nmf_cv_results.rds"))

set.seed(12)
cvnmf <- cross_validate_nmf(
    logcounts(sce),
    ranks=c(20, 30, 40, 50, 60),
    n_replicates = 1,
    tol = 1e-03,
    maxit = 100,
    verbose = 3,
    L1 = 0.1,
    L2 = 0,
    threads = 0,
    test_density = 0.2
)

pdf(here("plots", "18_PsychENCODE_NMF", "nmf_cv_results_lower_k.pdf"))
plot(cvnmf)
dev.off()

saveRDS(cvnmf, file = here("processed-data", "18_PsychENCODE_NMF", "nmf_cv_results_lower_k.rds"))
