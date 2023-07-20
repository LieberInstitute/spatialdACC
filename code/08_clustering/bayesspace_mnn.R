setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("BayesSpace")
    library("RColorBrewer")
    library("ggplot2")
    library("gridExtra")
    library("Polychrome")
})

load(here("processed-data", "07_batch_correction", "spe_mnn.Rdata"))
dim(spe)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

colData(spe)$row <- spe$array_row
colData(spe)$col <- spe$array_col

metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

message("Running spatialCluster()")
Sys.time()
set.seed(2)
spe <- spatialCluster(spe, use.dimred = "MNN", q = k,nrep=10000)
Sys.time()

bayesSpace_name <- paste0("bayesSpace_captureArea_", k)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
    spe,
    bayesSpace_name,
    cluster_dir = here::here("processed-data", "08_clustering", "BayesSpace", "preprocess_mnn", bayesSpace_name)
)

brains <- unique(spe$brnum)
cols <- Polychrome::palette36.colors(k)
names(cols) <- sort(unique(spe$spatial.cluster))
clustV <- bayesSpace_name

pdf(file = here::here("plots", "08_clustering", "BayesSpace", "preprocess_mnn", paste0(bayesSpace_name, ".pdf")), width = 21, height = 20)

for (i in seq_along(brains)){
    speb <- spe[, which(spe$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    print(length(samples))

    if (length(samples) == 1){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = clustV, colors = cols, point_size = 4, ... = paste0("_", brains[i]))
        grid.arrange(p1, nrow = 1)
    } else if (length(samples) == 2){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = clustV, colors = cols, point_size = 3, ... = paste0("_", brains[i]))
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = clustV, colors = cols, point_size = 3, ... = paste0("_", brains[i]))
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = clustV, colors = cols, point_size = 3, ... = paste0("_", brains[i]))
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = clustV, colors = cols, point_size = 3, ... = paste0("_", brains[i]))
        p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = clustV, colors = cols, point_size = 3, ... = paste0("_", brains[i]))
        grid.arrange(p1, p2, p3, nrow = 2)
    }
}
dev.off()
