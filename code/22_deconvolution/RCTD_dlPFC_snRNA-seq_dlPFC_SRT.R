setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library(here)
library(SpatialExperiment)
#BiocManager::install("spacexr")
#not GitHub version!
library(spacexr)

# set directories
plots_dir <- here("plots", "22_deconvolution")
processed_dir <- here("processed-data", "22_deconvolution")

# load
spe <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium")

spe_DLPFC_30.temp <- spe

# create spatial labels for DLPFC_30
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 3] <- "L2"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 8] <- "L4"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 7] <- "L6"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 5] <- "L3"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 6] <- "WM"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 4] <- "L5"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 2] <- "L1"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 1] <- "meninges"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 9] <- "WM"
# remove meninges
spe_DLPFC_30.temp <- spe_DLPFC_30.temp[ , which(spe_DLPFC_30.temp$BayesSpace_harmony_09 != "meninges")]

spe <- spe_DLPFC_30.temp
spe <- spe[!duplicated(rownames(spe)), ]

#remove duplications of colnames
colnames(spe) <- paste0(colnames(spe),colData(spe)$sample_id)

# load
sce_path_zip <- spatialLIBD::fetch_data("spatialDLPFC_snRNAseq")
sce_path <- unzip(sce_path_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
    file.path(tempdir(), "sce_DLPFC_annotated")
)

assay(sce, "logcounts") <- as(assay(sce, "logcounts"), "dgCMatrix")
sce

# make gene_ids row names to allow overlap
rownames(sce) <- rowData(sce)$gene_id

#remove / character 
library(stringr)
sce$cellType_layer <- str_replace(sce$cellType_layer, "/", "_")
sce$cellType_layer <- str_replace(sce$cellType_layer, "/", "_")

# ========== RCTD ==========
rctd_data <- createRctd(spe, sce, cell_type_col="cellType_layer")
res <- runRctd(rctd_data, max_cores=30, rctd_mode="multi", max_multi_types=5)

# save
saveRDS(res, here(processed_dir, "rctd_results_dlPFC_dlPFC.rds"))
