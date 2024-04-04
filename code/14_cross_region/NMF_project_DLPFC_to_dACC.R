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
spe_dACC <- spe_dACC[rowData(spe_dACC)$gene_name %in% rownames(spe_DLPFC),]

# drop rownames in loadings not in spe_dACC
loadings <- loadings[rownames(loadings) %in% rowData(spe_dACC)$gene_name,]

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

    pdf(file = here::here("plots", "13_NMF", "SpotPlots", paste0("NMF_", i, ".pdf")),
        width = 21, height = 20)

    for (j in seq_along(brains)){
        spe_dACCb <- spe_dACC.temp[, which(spe_dACC.temp$brnum == brains[j])]
        samples <- unique(spe_dACCb$sample_id)
        print(length(samples))

        if (length(samples) == 1){
            p1 <- vis_gene(spe_dACC =  spe_dACCb, sampleid = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 9, )
            grid.arrange(p1, nrow = 1)
        } else if (length(samples) == 2){
            p1 <- vis_gene(spe_dACC =  spe_dACCb, sampld = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p2 <- vis_gene(spe_dACC =  spe_dACCb, sampleid = samples[2], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            grid.arrange(p1, p2, nrow = 2)
        } else if (length(samples) == 3){
            p1 <- vis_gene(spe_dACC =  spe_dACCb, sampleid = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p2 <- vis_gene(spe_dACC =  spe_dACCb, sampleid = samples[2], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p3 <- vis_gene(spe_dACC =  spe_dACCb, sampleid = samples[3], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            grid.arrange(p1, p2, p3, nrow = 2)
        }
    }

    dev.off()
}

# find the number of zeroes in each column in colData(spe_dACC.temp)
zeroes <- sapply(colData(spe_dACC.temp)[, 1:100], function(x) sum(x == 0))

# list out columns in colData(spe_dACC.temp) with more than 77533 zeroes
names(zeroes[zeroes > 77533])

