setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SingleCellExperiment)
library(here)
library(scran)

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))
sce <- logNormCounts(sce)

marker.info <- scoreMarkers(sce, colData(sce)$cellType_azimuth)
marker.info

save(marker.info, file = here("processed-data", "snRNA-seq", "05_azimuth", "marker_info_azimuth.Rdata"))
