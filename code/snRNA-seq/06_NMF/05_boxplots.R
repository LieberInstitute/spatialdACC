setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

####code for correlating categorical technical variables with NMF patterns
library(here)
library(SpatialExperiment)
library(pheatmap)
library(reshape2)
library(RcppML)
library(ggplot2)
library(gridExtra)
library(spatialLIBD)
library(scater)

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))
x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))

sce.temp <- sce

# add each proj column to colData(sce)
# note that the rows of x@h already sum to 1 for each factor
for (i in 1:75){
    colData(sce.temp)[[paste0("NMF_",i)]] <- x@h[i,]
}

plot_list <- list()

for (i in 1:75){
    print(paste0("i=", i))

    p <- plotColData(sce.temp, x = "cellType_azimuth", y = paste0("NMF_", i)) +
        ggtitle(paste0("NMF ", i, " Layer Boxplots")) +
        facet_wrap(~ sce.temp$cellType_azimuth, scales = "free_x", nrow = 1)

    plot_list[[i]] <- p

}

for (i in seq(1, length(plot_list), by = 5)) {
    print(i)

    pdf(file = here::here("plots", "snRNA-seq", "06_NMF", paste0("NMF_boxplots_", i, "-", i+4, ".pdf")),
        width = 21, height = 20)

    grid.arrange(
        grobs = plot_list[i:min(i+4, length(plot_list))],
        ncol = 1,
        nrow = 5
    )

    dev.off()

}


plot_list <- list()

for (i in c(38,61)){
    print(paste0("i=", i))

    p <- plotColData(sce.temp, x = "cellType_azimuth", y = paste0("NMF_", i)) +
        ggtitle(paste0("NMF ", i, " dACC snRNA-seq Layer Boxplots")) +
        facet_wrap(~ sce.temp$cellType_azimuth, scales = "free_x", nrow = 1) +
        ylim(c(0,0.0035))

    plot_list[[i]] <- p

}

pdf(file = here::here("plots", "snRNA-seq", "06_NMF", "NMF_boxplots_dACC_single_nucleus_38_61.pdf"),
    width = 13, height = 10)

grid.arrange(
    grobs = plot_list[c(38,61)],
    ncol = 1,
    nrow = 2
)

dev.off()
