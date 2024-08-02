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

# get NMF results from single nucleus data
x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))

# load uncorrected and lognormalized spatial data
load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))

# load Single Nucleus object
# this file says "harmony", but that is only in the reduced dims, the counts/logcounts are not batch corrected
load(file = here("processed-data", "snRNA-seq", "03_batch_correction", "sce_harmony.Rdata"))

# extract patterns
patterns <- t(x@h)
colnames(patterns) <- paste("NMF", 1:100, sep = "_")

loadings <- x@w
rownames(loadings) <- rownames(sce)

# ====== project loadings to spatial data =======
i <- intersect(rowData(spe)$gene_name,rownames(loadings))
loadings <- loadings[rownames(loadings) %in% i,]
spe <- spe[rowData(spe)$gene_name %in% i,]
loadings <- loadings[match(rowData(spe)$gene_name,rownames(loadings)),]

logcounts <- logcounts(spe)
#loadings <- as(loadings, "dgCMatrix")

proj <- project(w=loadings, data=logcounts)
proj <- t(proj)
colnames(proj) <- paste("NMF", 1:100, sep = "_")

# add to reducedDims
reducedDim(spe, "NMF_proj") <- proj

# save spe
save(spe, file = here("processed-data", "13_NMF", "spe_NMF.Rdata"))

spe.temp <- spe

# add each proj column to colData(spe)
for (i in 1:100){
    colData(spe.temp)[[paste0("NMF_",i)]] <- reducedDims(spe.temp)$NMF_proj[,i]
}

brains <- unique(spe.temp$brnum)

for (i in 1:100){
    print(paste0("i=", i))

    pdf(file = here::here("plots", "13_NMF", "SpotPlots", paste0("NMF_", i, ".pdf")),
        width = 21, height = 20)

    for (j in seq_along(brains)){
        speb <- spe.temp[, which(spe.temp$brnum == brains[j])]
        samples <- unique(speb$sample_id)
        print(length(samples))

        if (length(samples) == 1){
            p1 <- vis_gene(spe =  speb, sampleid = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 9, )
            grid.arrange(p1, nrow = 1)
        } else if (length(samples) == 2){
            p1 <- vis_gene(spe =  speb, sampleid = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p2 <- vis_gene(spe =  speb, sampleid = samples[2], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            grid.arrange(p1, p2, nrow = 2)
        } else if (length(samples) == 3){
            p1 <- vis_gene(spe =  speb, sampleid = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p2 <- vis_gene(spe =  speb, sampleid = samples[2], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p3 <- vis_gene(spe =  speb, sampleid = samples[3], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            grid.arrange(p1, p2, p3, nrow = 2)
        }
    }

    dev.off()
}

# add PRECAST clusters to spe
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9.Rdata"))

spe$PRECAST_cluster <- unfactor(spe$PRECAST_cluster)
spe$PRECAST_cluster[spe$PRECAST_cluster == 3] <- "WM1"
spe$PRECAST_cluster[spe$PRECAST_cluster == 8] <- "WM2"
spe$PRECAST_cluster[spe$PRECAST_cluster == 7] <- "WM-CC"
spe$PRECAST_cluster[spe$PRECAST_cluster == 5] <- "L6b"
spe$PRECAST_cluster[spe$PRECAST_cluster == 6] <- "L6a"
spe$PRECAST_cluster[spe$PRECAST_cluster == 4] <- "L5"
spe$PRECAST_cluster[spe$PRECAST_cluster == 2] <- "L3"
spe$PRECAST_cluster[spe$PRECAST_cluster == 1] <- "L2"
spe$PRECAST_cluster[spe$PRECAST_cluster == 9] <- "L1"

spe.temp$PRECAST_cluster <- spe$PRECAST_cluster

plot_list <- list()

for (i in 1:100){
    print(paste0("i=", i))

    p <- plotColData(spe.temp, x = "PRECAST_cluster", y = paste0("NMF_", i)) +
        ggtitle(paste0("NMF ", i, " Layer Boxplots")) +
        facet_wrap(~ spe.temp$PRECAST_cluster, scales = "free_x", nrow = 1)

    plot_list[[i]] <- p

}

for (i in seq(1, length(plot_list), by = 5)) {
    print(i)

    pdf(file = here::here("plots", "13_NMF", paste0("NMF_boxplots_", i, "-", i+4, ".pdf")),
        width = 21, height = 20)

    grid.arrange(
        grobs = plot_list[i:min(i+4, length(plot_list))],
        ncol = 1,
        nrow = 5
    )

    dev.off()

}
