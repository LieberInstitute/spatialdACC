setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("scry"))
suppressPackageStartupMessages(library("BiocSingular"))
suppressPackageStartupMessages(library("schex"))
suppressPackageStartupMessages(library("PCAtools"))
suppressPackageStartupMessages(library("glmpca"))

load(here("processed-data", "snRNA-seq", "01_QC", "sce_doublet.rda"))

sce$brain <- as.factor(sce$brain)

dim(sce)
# 34866 42174

#remove doublets
sce$high_doublet <- sce$scDblFinder.score > 0.99
sce <- sce[, !colData(sce)$high_doublet]

# poisson deviance feature selection, then GLM PCA
set.seed(8)
sce <- devianceFeatureSelection(sce, assay = "counts", fam = "poisson", sorted = T, batch = sce$brain)

pdf(here("plots", "snRNA-seq", "02_preprocessing", "poisson_deviance_doublets_removed.pdf"))
plot(sort(rowData(sce)$poisson_deviance, decreasing = T),
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

#glm pca
hdg <- rownames(counts(sce))[1:2000]

set.seed(9)
message("running nullResiduals - ", Sys.time())
res <- sce[rownames(counts(sce)) %in% hdg,]
res <- nullResiduals(res,
                     fam = "poisson",
                     type = "pearson",
                     assay = "counts"   #, batch = res$brain giving errors
)

set.seed(10)
message("running PCA - ", Sys.time())

res <- scater::runPCA(res, ncomponents = 50,
                      exprs_values='poisson_pearson_residuals',
                      scale = TRUE, name = "pp-GLM-PCA",
                      BSPARAM = BiocSingular::RandomParam())

#pca plots
pdf(here("plots", "snRNA-seq", "02_preprocessing", "sel_poisson_pearson_GLM_PCA_brain_doublets_removed.pdf"))
hex <- make_hexbin(res, nbins = 100, dimension_reduction = "pp-GLM-PCA", use_dims = c(1, 2))
label_df <- make_hexbin_label(hex, col = "brain")
plot_hexbin_meta(hex, col = "brain", action = "majority", xlab = "PC1", ylab = "PC2") + ggtitle("Brain") + theme(legend.position = "right")

res$discard_auto <- as.logical(res$discard_auto)
res$discard_auto <- as.factor(res$discard_auto)
hex <- make_hexbin(res, nbins = 100, dimension_reduction = "pp-GLM-PCA", use_dims = c(1, 2))
label_df <- make_hexbin_label(hex, col = "discard_auto")
plot_hexbin_meta(hex, col = "discard_auto", action = "majority", xlab = "PC1", ylab = "PC2") + ggtitle("Low Quality") + theme(legend.position = "right")

hex <- make_hexbin(res, nbins = 100, dimension_reduction = "pp-GLM-PCA", use_dims = c(1, 2))
label_df <- make_hexbin_label(hex, col = "scDblFinder.class")
plot_hexbin_meta(hex, col = "scDblFinder.class", action = "majority", xlab = "PC1", ylab = "PC2") + ggtitle("Doublet") + theme(legend.position = "right")

dev.off()

# make elbow plot to determine PCs to use
percent.var <- attr(reducedDim(res, "pp-GLM-PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow

pdf(here("plots", "snRNA-seq", "02_preprocessing", "sel_poisson_pearson_var_explained_doublets_removed.pdf"))
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()

#run umap
reducedDim(sce,'pp-GLM-PCA') <- reducedDim(res,'pp-GLM-PCA')

message("Running runUMAP()")
Sys.time()
set.seed(11)
sce <- runUMAP(sce, dimred = "pp-GLM-PCA", name="UMAP-GLM-PCA")
Sys.time()

#explore UMAP results
pdf(file = here::here("plots", "snRNA-seq", "02_preprocessing", "sel_poisson_pearson_UMAP_doublets_removed.pdf"))

hex <- make_hexbin(sce, nbins = 100, dimension_reduction = "UMAP-GLM-PCA", use_dims = c(1, 2))
label_df <- make_hexbin_label(hex, col = "brain")
plot_hexbin_meta(hex, col = "brain", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("Brain") + theme(legend.position = "right")

sce$discard_auto <- as.logical(sce$discard_auto)
sce$discard_auto <- as.factor(sce$discard_auto)
hex <- make_hexbin(sce, nbins = 100, dimension_reduction = "UMAP-GLM-PCA", use_dims = c(1, 2))
label_df <- make_hexbin_label(hex, col = "discard_auto")
plot_hexbin_meta(hex, col = "discard_auto", action = "majority", xlab = "PC1", ylab = "PC2") + ggtitle("Low Quality") + theme(legend.position = "right")

hex <- make_hexbin(sce, nbins = 100, dimension_reduction = "UMAP-GLM-PCA", use_dims = c(1, 2))
label_df <- make_hexbin_label(hex, col = "scDblFinder.class")
plot_hexbin_meta(hex, col = "scDblFinder.class", action = "majority", xlab = "PC1", ylab = "PC2") + ggtitle("Doublet") + theme(legend.position = "right")

dev.off()

save(sce, file = here::here("processed-data", "snRNA-seq", "02_preprocessing", "sce_dimred_doublets_removed.Rdata"))
