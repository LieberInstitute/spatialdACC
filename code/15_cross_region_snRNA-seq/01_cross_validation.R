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
sce_path_zip <- fetch_data("spatialDLPFC_snRNAseq")
sce_path <- unzip(sce_path_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
    file.path(tempdir(), "sce_DLPFC_annotated")
)

assay(sce, "logcounts") <- as(assay(sce, "logcounts"), "dgCMatrix")

#cvnmf <- cross_validate_nmf(
#    logcounts(sce),
#    ranks=c(25,50,75,100,125),
#    n_replicates = 2,
#    tol = 1e-03,
#    maxit = 100,
#    verbose = 3,
#    L1 = 0.1,
#    L2 = 0,
#    threads = 0,
#    test_density = 0.2
#)

#saveRDS(cvnmf, file = here("processed-data", "15_cross_region_snRNA-seq", "DLPFC_nmf_cv_results.RDS"))

#pdf(here("plots", "15_cross_region_snRNA-seq", "DLPFC_nmf_cv_results.pdf"))
#plot(cvnmf)
#dev.off()

set.seed(78)
cvnmf <- cross_validate_nmf(
    logcounts(sce),
    ranks=c(5, 10, 15, 20, 25, 30, 50, 75, 100),
    n_replicates = 2,
    tol = 1e-03,
    maxit = 100,
    verbose = 3,
    L1 = 0.1,
    L2 = 0,
    threads = 0,
    test_density = 0.2
)

pdf(here("plots", "15_cross_region_snRNA-seq", "DLPFC_nmf_cv_results_lower_k.pdf"))
plot(cvnmf)
dev.off()

saveRDS(cvnmf, file = here("processed-data", "15_cross_region_snRNA-seq", "DLPFC_nmf_cv_results_lower_k.RDS"))
