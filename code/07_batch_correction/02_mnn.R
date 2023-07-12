setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("schex"))
suppressPackageStartupMessages(library("batchelor"))
suppressPackageStartupMessages(library("ggplot2"))

load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))

reducedDimNames(spe)
# "10x_pca"      "10x_tsne"     "10x_umap"     "pp-GLM-PCA"   "UMAP-GLM-PCA"

mnn <- reducedMNN(reducedDim(spe,'pp-GLM-PCA'), batch=spe$brnum, k=50)

reducedDim(spe,'MNN') <- mnn$corrected
set.seed(125)
message("running UMAP - ", Sys.time())
spe <- runUMAP(spe, dimred = "MNN", name = "UMAP-MNN")

#explore UMAP results
pdf(file = here::here("plots", "07_batch_correction", "sel_poisson_pearson_MNN_UMAP.pdf"))

hex <- make_hexbin(spe, nbins = 100, dimension_reduction = "UMAP-MNN", use_dims = c(1, 2))
label_df <- make_hexbin_label(hex, col = "brnum")
plot_hexbin_meta(hex, col = "brnum", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("Brains") + theme(legend.position = "right")

label_df <- make_hexbin_label(hex, col = "sample_id")
plot_hexbin_meta(hex, col = "sample_id", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("Capture area") + theme(legend.position = "right")

dev.off()

save(spe, file = here::here("processed-data", "07_batch_correction", "spe_mnn.Rdata"))

