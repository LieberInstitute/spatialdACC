setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("scry"))
suppressPackageStartupMessages(library("BiocSingular"))

load(here("processed-data", "06_preprocessing", "spe_norm.Rdata"))

#option 1: poisson deviance feature selection, then GLM PCA
set.seed(8)
spe <- devianceFeatureSelection(spe, assay = "counts", fam = "poisson", sorted = T)

pdf(here("plots", "06_preprocessing", "poisson_deviance.pdf"))
plot(sort(rowData(spe)$poisson_deviance, decreasing = T),
     type = "l", xlab = "ranked genes",
     ylim=c(0,1000000),
     ylab = "poisson deviance", main = "Feature Selection with Deviance"
)
abline(v = 1000, lty = 2, col = "purple")
abline(v = 2000, lty = 2, col = "red")
abline(v = 2500, lty = 2, col = "pink")
abline(v = 3000, lty = 2, col = "blue")
abline(v = 4000, lty = 2, col = "green")
abline(v = 5000, lty = 2, col = "black")
dev.off()

hdg <- rownames(counts(spe))[1:1000]

set.seed(9)
message("running nullResiduals - ", Sys.time())
res <- spe[rownames(counts(spe)) %in% hdg,]
res <- nullResiduals(res,
                     fam = "poisson",
                     type = "pearson",
                     assay='counts'
)

set.seed(10)
message("running PCA - ", Sys.time())

res <- scater::runPCA(res, ncomponents = 50,
                      exprs_values='poisson_pearson_residuals',
                      scale = TRUE, name = "pp-GLM-PCA",
                      BSPARAM = BiocSingular::RandomParam())

pdf(here("plots", "06_preprocessing", "sel_poisson_pearson_GLM_PCA_brnum.pdf"))
plotReducedDim(res, dimred = "pp-GLM-PCA", colour_by = "brnum")
dev.off()

pdf(here("plots", "06_preprocessing", "sel_poisson_pearson_GLM_PCA_sample_id.pdf"))
plotReducedDim(res, dimred = "pp-GLM-PCA", colour_by = "sample_id")
dev.off()

#option 2: binomial deviance feature selection, then GLM PCA
set.seed(8)
spe <- devianceFeatureSelection(spe, assay = "counts", fam = "binomial", sorted = T)

pdf(here("plots", "06_preprocessing", "binomial_deviance.pdf"))
plot(sort(rowData(spe)$binomial_deviance, decreasing = T),
     type = "l", xlab = "ranked genes",
     ylim=c(0,1000000),
     ylab = "binomial deviance", main = "Feature Selection with Deviance"
)
abline(v = 1000, lty = 2, col = "purple")
abline(v = 2000, lty = 2, col = "red")
abline(v = 2500, lty = 2, col = "pink")
abline(v = 3000, lty = 2, col = "blue")
abline(v = 4000, lty = 2, col = "green")
abline(v = 5000, lty = 2, col = "black")
dev.off()

hdg <- rownames(counts(spe))[1:1000]

set.seed(9)
message("running nullResiduals - ", Sys.time())
res <- spe[rownames(counts(spe)) %in% hdg,]
res <- nullResiduals(res,
                     fam = "binomial",
                     type = "deviance",
                     assay='counts'
)

set.seed(10)
res <- scater::runPCA(res, ncomponents = 50,
                      exprs_values='binomial_deviance_residuals',
                      scale = TRUE, name = "bd-GLM-PCA",
                      BSPARAM = BiocSingular::RandomParam())

pdf(here("plots", "06_preprocessing", "binomial_deviance_GLM_PCA_brnum.pdf"))
plotReducedDim(res, dimred = "bd-GLM-PCA", colour_by = "brnum")
dev.off()

pdf(here("plots", "06_preprocessing", "binomial_deviance_GLM_PCA_sample_id.pdf"))
plotReducedDim(res, dimred = "bd-GLM-PCA", colour_by = "sample_id")
dev.off()
