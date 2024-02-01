library("SingleCellExperiment")
library("harmony")
library("scater")
library("here")
library("sessioninfo")

## load sce
load(file = here::here("processed-data", "snRNA-seq", "02_preprocessing", "sce_dimred_doublets_removed.Rdata"))

## Run harmony
## needs PCA
reducedDim(sce, "PCA") <- reducedDim(sce, "pp-GLM-PCA")

message("running Harmony - ", Sys.time())
sce <- RunHarmony(sce, group.by.vars = "Sample", verbose = TRUE)

## Remove redundant PCA
reducedDim(sce, "PCA") <- NULL

#### TSNE & UMAP ####
set.seed(602)

message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "HARMONY", name="UMAP-HARMONY")

pdf(file = here::here("plots", "snRNA-seq", "03_batch_correction", "HARMONY_UMAP.pdf"))
plotReducedDim(sce, dimred="UMAP-HARMONY", colour_by="Sample", point_size = 0.5)
dev.off()

save(sce, file = here("processed-data", "snRNA-seq", "03_batch_correction", paste0("sce_harmony.Rdata")))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
