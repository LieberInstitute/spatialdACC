setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
# Load libraries
library(RcppML)
library(SingleCellExperiment)
library(here)
library(scuttle)
library(Matrix)
library(singlet)

# want to use logcounts for NMF
sce <- readRDS(file = here("processed-data", "18_PsychENCODE_NMF", "pseudobulk", "pseudobulk_combined.rds"))

set.seed(1234)
res <- RunNMF(
    sce,
    k=c(25,50,75,100,125,150),
    assay = "logcounts",
    reps = 3,
    tol = 1e-03,
    maxit = 100,
    verbose = 3,
    L1 = 0.1,
    L2 = 0,
    threads = 0,
    test.set.density = 0.2
)

saveRDS(res, file = here("processed-data", "18_PsychENCODE_NMF", "nmf_cv_results.rds"))

pdf(here("plots", "18_PsychENCODE_NMF", "nmf_cv_results.pdf"))
RankPlot(res)
dev.off()

