setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("PRECAST")
    library("tictoc")
    library("dplyr")
    library("purrr")
    library("tidyverse")
    library("spatialLIBD")
    library("gridExtra")
    library("ggspavis")
    library("patchwork")
})

load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))
colnames(spe) <- spe$key

for (i in c(5:20)) {
    print(i)
    precast_name <- paste0("nnSVG_PRECAST_captureArea_", i)
    spe <- cluster_import(
        spe,
        cluster_dir = here::here("processed-data", "08_clustering", "PRECAST", precast_name),
        prefix = paste0("clust", i, "_")
    )
}

# vis the PRECAST clusters from sample V12N28-334_C1
plot_list <- list()

for (i in c(5:20)) {
    print(i)

    p <- vis_clus(spe = spe, sampleid = "V12N28-334_C1",
                  clustervar = paste0("clust",i,"_PRECAST_cluster"), point_size = 2,
                  spatial=F) +
        guides(fill = guide_legend(override.aes = list(size = 5))) +
        theme(legend.text = element_text(size = 10)) +
        ggtitle("")

    plot_list[[i]] <- p

}


png(here("plots", "08_clustering", "nnSVG_PRECAST_5_to_20.png"), width = 25, height = 25, units = "in", res=300)

wrap_plots(plot_list[[5]],
           plot_list[[6]],
           plot_list[[7]],
           plot_list[[8]],
           plot_list[[9]],
           plot_list[[10]],
           plot_list[[11]],
           plot_list[[12]],
           plot_list[[13]],
           plot_list[[14]],
           plot_list[[15]],
           plot_list[[16]],
           plot_list[[17]],
           plot_list[[18]],
           plot_list[[19]],
           plot_list[[20]],
           nrow = 4) +
    plot_annotation(title="nnSVG-Guided PRECAST Clusters k=5-20 (V12N28-334_C1)",theme=theme(plot.title = element_text(size = 40)))

dev.off()
