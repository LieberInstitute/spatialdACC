setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library(SpatialExperiment)
library(scran)
library(escheR)
library(ggplot2)
library(here)
library(Seurat)
library(spatialLIBD)
library(nmfLabelTransfer)

load(here("processed-data", "VSPG", "01_QC", "spe_QC.Rdata"), verbose = TRUE)

# remove spots with low sum and low detected features
spe <- spe[,!colData(spe)$low_sum_br & !colData(spe)$low_detected_br]
dim(spe)
# [1] 28096 17278

spe_target <- spe
spe_target <- logNormCounts(spe_target)

# load spe for k=9 without WM-CC
load(
    file = here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata")
)

spe_source <- spe
dim(spe)

source_seu <- as.Seurat(spe_source)
source_seu <- NormalizeData(source_seu)
source_seu <- FindVariableFeatures(source_seu)
source_seu <- ScaleData(source_seu)
source_seu <- RunPCA(source_seu)

# Convert target datasets to Seurat and preprocess
colnames(spe_target) <- make.unique(colnames(spe_target))
target_seu <- as.Seurat(spe_target)
target_seu <- NormalizeData(target_seu)

print(source_seu)
source_seu <- AddMetaData(source_seu, metadata = as.factor(spe_source$layer),
                          col.name="layer")

query_seu <- target_seu

# Find anchors between ref and query
anchors <- FindTransferAnchors(reference = source_seu,
                               query =  query_seu, dims = 1:30,
                               normalization.method = "LogNormalize",
                               reference.reduction = "pca")
# Predict annotations
seu_predictions <- TransferData(anchorset = anchors,
                                refdata = source_seu$layer,
                                dims = 1:30)
print(class(seu_predictions))
print(head(seu_predictions))
query_seu <- AddMetaData(query_seu, metadata = seu_predictions)
print(class(query_seu[["predicted.id"]]))
preds <- query_seu[["predicted.id"]]
seu_predictions <- preds[,"predicted.id"]
seu_predictions <- as.factor(seu_predictions)

spe_target$seu_predictions <- seu_predictions

# Save the annotated object
saveRDS(spe_target, here("processed-data", "VSPG", "02_label_transfer", "seurat_target_with_preds.rds"))

brains <- unique(spe_target$brnum)
plot_list <- c()
# use this code but make a separate plot for each brain
for (i in seq_along(brains)){

    speb <- spe_target[, which(spe_target$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    print(length(samples))

    p <- make_escheR(speb, y_reverse=FALSE) %>%
        add_fill("seu_predictions", point_size = 1)+
        scale_fill_manual(values = c(
            "L2" = "#E41A1C",   # Bright red
            "L3" = "#377EB8",   # Strong blue
            "L5" = "#4DAF4A",   # Vivid green
            "L6a" = "#984EA3",  # Purple
            "L6b" = "#FF7F00",  # Orange
            "WM" = "#F781BF",# Pink
            "L1" = "#00CED1"    # Dark turquoise
        ))

    plot_list[[i]] <- p
}

pdf(here("plots", "VSPG", "02_label_transfer", "seurat_target_with_preds.pdf"))
do.call(gridExtra::grid.arrange, c(plot_list, ncol=2))
dev.off()
