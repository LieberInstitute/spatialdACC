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
library("harmony")
library("patchwork")

load(here("processed-data", "snRNA-seq", "01_QC", "sce_doublet.rda"))

sce$brain <- as.factor(sce$brain)

dim(sce)
# 34866 37032

#remove doublets
sce <- sce[, colData(sce)$scDblFinder.class == "singlet"]
dim(sce)
# [1] 34866 35161

# poisson deviance feature selection, then GLM PCA
set.seed(8)
sce <- devianceFeatureSelection(sce, assay = "counts", fam = "poisson", sorted = T, batch = sce$brain)

load(file = here("processed-data", "snRNA-seq", "03_batch_correction", paste0("sce_harmony.Rdata")))

poisson_dev <- sort(rowData(sce)$poisson_deviance, decreasing = TRUE)
df <- data.frame(
    rank = seq_along(poisson_dev),
    deviance = poisson_dev
)

p1 <- ggplot(df, aes(x = rank, y = deviance)) +
    geom_line() +
    scale_y_continuous(limits = c(0, 1e6)) +
    geom_vline(xintercept = c(1000, 2000, 2500, 3000, 4000, 5000),
               linetype = "dashed",
               color = c("purple", "red", "pink", "blue", "green", "black")) +
    labs(
        x = "ranked genes",
        y = "poisson deviance",
        title = "Feature Selection with Deviance"
    ) +
    theme_minimal()

p2 <- plotReducedDim(sce, dimred="UMAP-HARMONY", colour_by="Sample", point_size = 0.2)

png(here("plots", "snRNA-seq", "03_batch_correction", "preprocessing_supp.png"), height = 4, width=8, unit="in", res=300)

wrap_plots(p1,p2,nrow=1) + plot_annotation(tag_levels = 'A')
dev.off()
