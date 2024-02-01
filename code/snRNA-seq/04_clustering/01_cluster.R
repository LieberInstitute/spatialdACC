library("SingleCellExperiment")
library("jaffelab")
library("scater")
library("scran")
library("here")
library("sessioninfo")

#load sce
load(file = here("processed-data", "snRNA-seq", "03_batch_correction", "sce_harmony.Rdata"))

message("running buildSNNGraph - ", Sys.time())
snn.gr <- buildSNNGraph(sce, k = 20, use.dimred = "HARMONY")

message("running walktrap - ", Sys.time())
clusters <- igraph::cluster_walktrap(snn.gr)$membership

table(clusters)
message("saving data - ", Sys.time())
save(clusters, file = here("processed-data", "snRNA-seq", "04_clustering", "HARMONY_clusters.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
