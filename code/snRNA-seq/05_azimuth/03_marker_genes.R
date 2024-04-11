setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SingleCellExperiment)
library(here)
library(scran)

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))
sce <- logNormCounts(sce)

# replace spaces with underscores
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L2/3 IT", "L2_3_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5 ET", "L5_ET")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5 IT", "L5_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5/6 NP", "L5_6_NP")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 CT", "L6_CT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 IT", "L6_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 IT Car3", "L6_IT_Car3")

marker.info <- scoreMarkers(sce, colData(sce)$cellType_azimuth)
marker.info

save(marker.info, file = here("processed-data", "snRNA-seq", "05_azimuth", "marker_info_azimuth.Rdata"))

# get top 20 markers for each cluster
top20_markers <- lapply(marker.info, function(x) {
  x <- x[order(x$mean.AUC, decreasing = TRUE),]
  x <- x[1:20,]
  return(x)
})

# write to csv file
for (i in 1:length(top20_markers)) {
  write.csv(top20_markers[[i]], file = here("processed-data", "snRNA-seq", "05_azimuth", paste0("top20_markers_", names(top20_markers)[i], ".csv")))
}
