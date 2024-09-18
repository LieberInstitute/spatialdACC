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
