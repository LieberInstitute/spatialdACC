setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("scry"))
suppressPackageStartupMessages(library("BiocSingular"))

load(here("processed-data", "06_preprocessing", "spe_norm.Rdata"))

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
                      scale = TRUE, name = "GLM-PCA",
                      BSPARAM = BiocSingular::RandomParam())

pdf(here("plots", "06_preprocessing", "GLM_PCA_brnum.pdf"))
plotReducedDim(res, dimred = "GLM-PCA", colour_by = "brnum")
dev.off()

pdf(here("plots", "06_preprocessing", "GLM_PCA_sample_id.pdf"))
plotReducedDim(res, dimred = "GLM-PCA", colour_by = "sample_id")
dev.off()
