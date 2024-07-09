setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(scater)
library(RcppML)
library(ggspavis)
library(here)
library(scRNAseq)
library(Matrix)
library(scran)
library(scuttle)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(igraph)
library(bluster)
library(patchwork)
library(cowplot)
library(projectR)
library(spatialLIBD)
library(gridExtra)

# get NMF results from DLPFC data
x <- readRDS(file = here("processed-data", "14_cross_region", "DLPFC_nmf_results.RDS"))

# load uncorrected and lognormalized dACC spatial data
load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))
spe_dACC <- spe

# load DLPFC spatial object
spe_DLPFC <- spatialLIBD::fetch_data(type = "spe")

# extract patterns
patterns <- t(x$h)
colnames(patterns) <- paste("NMF", 1:100, sep = "_")

loadings <- x$w
rownames(loadings) <- rownames(spe_DLPFC)

# ====== project loadings to spatial data =======
# drop any gene_names in spe_dACC not in speLPFC
spe_dACC <- spe_dACC[rowData(spe_dACC)$gene_id %in% rownames(spe_DLPFC),]

# drop rownames in loadings not in spe_dACC
loadings <- loadings[rownames(loadings) %in% rowData(spe_dACC)$gene_id,]
loadings <- loadings[match(rownames(spe_dACC), rownames(loadings)), ]

logcounts <- logcounts(spe_dACC)

proj <- project(loadings, as.matrix(logcounts))
proj <- t(proj)
colnames(proj) <- paste("NMF", 1:100, sep = "_")

# add to reducedDims
reducedDim(spe_dACC, "NMF_proj") <- proj

# save spe_dACC
save(spe_dACC, file = here("processed-data", "14_cross_region", "spe_dACC_NMF.Rdata"))

spe_dACC.temp <- spe_dACC

# add each proj column to colData(spe_dACC.temp)
for (i in 1:100){
    colData(spe_dACC.temp)[[paste0("NMF_",i)]] <- reducedDims(spe_dACC.temp)$NMF_proj[,i]
}

brains <- unique(spe_dACC.temp$brnum)

for (i in 1:100){
    print(paste0("i=", i))

    pdf(file = here::here("plots", "14_cross_region", "SpotPlots", paste0("NMF_", i, ".pdf")),
        width = 21, height = 20)

    for (j in seq_along(brains)){
        spe_dACCb <- spe_dACC.temp[, which(spe_dACC.temp$brnum == brains[j])]
        samples <- unique(spe_dACCb$sample_id)
        print(length(samples))

        if (length(samples) == 1){
            p1 <- vis_gene(spe =  spe_dACCb, sampleid = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 9, )
            grid.arrange(p1, nrow = 1)
        } else if (length(samples) == 2){
            p1 <- vis_gene(spe =  spe_dACCb, sampld = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p2 <- vis_gene(spe =  spe_dACCb, sampleid = samples[2], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            grid.arrange(p1, p2, nrow = 2)
        } else if (length(samples) == 3){
            p1 <- vis_gene(spe =  spe_dACCb, sampleid = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p2 <- vis_gene(spe =  spe_dACCb, sampleid = samples[2], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p3 <- vis_gene(spe =  spe_dACCb, sampleid = samples[3], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            grid.arrange(p1, p2, p3, nrow = 2)
        }
    }

    dev.off()
}

# add PRECAST clusters to spe
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9.Rdata"))

spe$PRECAST_cluster <- unfactor(spe$PRECAST_cluster)
spe$PRECAST_cluster[spe$PRECAST_cluster == 3] <- "WM"
spe$PRECAST_cluster[spe$PRECAST_cluster == 8] <- "WM"
spe$PRECAST_cluster[spe$PRECAST_cluster == 7] <- "WM"
spe$PRECAST_cluster[spe$PRECAST_cluster == 5] <- "L6"
spe$PRECAST_cluster[spe$PRECAST_cluster == 6] <- "L6"
spe$PRECAST_cluster[spe$PRECAST_cluster == 4] <- "L5"
spe$PRECAST_cluster[spe$PRECAST_cluster == 2] <- "L3"
spe$PRECAST_cluster[spe$PRECAST_cluster == 1] <- "L2"
spe$PRECAST_cluster[spe$PRECAST_cluster == 9] <- "L1"

spe_dACC.temp$PRECAST_cluster <- spe$PRECAST_cluster

plot_list <- list()

for (i in 1:100){
    print(paste0("i=", i))

    p <- plotColData(spe_dACC.temp, x = "PRECAST_cluster", y = paste0("NMF_", i)) +
        ggtitle(paste0("NMF ", i, " Layer Boxplots")) +
        facet_wrap(~ spe_dACC.temp$PRECAST_cluster, scales = "free_x", nrow = 1) +
        labs(x = "Layer", y = paste0("NMF_", i)) +
        theme_bw()

    plot_list[[i]] <- p

}

for (i in seq(1, length(plot_list), by = 5)) {
    print(i)

    pdf(file = here::here("plots", "14_cross_region", paste0("NMF_boxplots_dACC_", i, "-", i+4, ".pdf")),
        width = 21, height = 20)

    grid.arrange(
        grobs = plot_list[i:min(i+4, length(plot_list))],
        ncol = 1,
        nrow = 5
    )

    dev.off()

}

# find the number of zeroes in each column in colData(spe_dACC.temp)
zeroes <- sapply(colData(spe_dACC.temp)[, 1:100], function(x) sum(x == 0))

# list out columns in colData(spe_dACC.temp) with more than 77533 zeroes
names(zeroes[zeroes > 77533])

# create a loop to run kruskal.test and rank_epsilon_squared for all NMFs
library(effectsize)

kruskal_res_list <- list()
epsilon_squared_res_list <- list()

for (i in 1:100){

    NMF_i <- colData(spe_dACC.temp)[, paste0("NMF_", i)]
    PRECAST_cluster <- spe_dACC.temp$PRECAST_cluster
    dat <- data.frame(NMF_i, PRECAST_cluster)

    # Check if all values of NMF_i are zero
    if (all(NMF_i == 0)) {
        kruskal_res_list[[i]] <- NA
        epsilon_squared_res_list[[i]] <- NA
        next  # Skip to the next iteration
    }

    kruskal_res <- kruskal.test(NMF_i ~ PRECAST_cluster, data = dat)
    epsilon_squared_res <- rank_epsilon_squared(kruskal_res, data = dat)

    kruskal_res_list[[i]] <- kruskal_res
    epsilon_squared_res_list[[i]] <- epsilon_squared_res

    if (epsilon_squared_res$rank_epsilon_squared > 0.16) {
        print(paste0("NMF_", i))
    }
}

# create csv with epsilon squared results and kruskal results

#extract only epsilon_squared_res$rank_epsilon_squared from epsilon_squared_res_list
epsilon_squared_values <- sapply(epsilon_squared_res_list, function(x) {
    if (is.null(x) || all(is.na(x))) {
        return(NA)
    } else {
        return(x$rank_epsilon_squared)
    }
})

# extract only p.value and statistic from kruskal_res_list
# Extract p.value and statistic from kruskal_res_list
kruskal_p_values <- sapply(kruskal_res_list, function(x) {
    if (is.null(x) || all(is.na(x))) {
        return(NA)
    } else {
        return(x$p.value)
    }
})

kruskal_statistics <- sapply(kruskal_res_list, function(x) {
    if (is.null(x) || all(is.na(x))) {
        return(NA)
    } else {
        return(x$statistic)
    }
})

# Combine the extracted values into a data frame
kruskal_results <- data.frame(
    NMF = paste0("NMF_", 1:100),
    p.value = kruskal_p_values,
    statistic = kruskal_statistics,
    epsilon_squared_values
)

write.csv(kruskal_results, here::here("processed-data", "14_cross_region", "kruskal_NMF_results_dACC.csv"))



