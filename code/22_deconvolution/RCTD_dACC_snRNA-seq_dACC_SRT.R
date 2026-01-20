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
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata"))
spe <- spe[!duplicated(rownames(spe)), ]

# load
load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))
sce

# make gene_ids row names to allow overlap
rownames(sce) <- rowData(sce)$gene_id


# ========== RCTD ==========
rctd_data <- createRctd(spe, sce, cell_type_col="cellType_azimuth")
res <- runRctd(rctd_data, max_cores=30, rctd_mode="multi", max_multi_types=5)

# save
saveRDS(res, here(processed_dir, "rctd_results.rds"))
