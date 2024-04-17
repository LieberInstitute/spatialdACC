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

# load new DLPFC spatial object
spe_DLPFC_30 <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium")

# load old DLPFC spatial object
spe_DLPFC_12 <- spatialLIBD::fetch_data(type = "spe")

# extract patterns
patterns <- t(x$h)
colnames(patterns) <- paste("NMF", 1:100, sep = "_")

loadings <- x$w
rownames(loadings) <- rownames(spe_DLPFC_12)

# ====== project loadings to spatial data =======
# drop any gene_names in spe_DLPFC_30 not in speLPFC
spe_DLPFC_30 <- spe_DLPFC_30[rowData(spe_DLPFC_30)$gene_id %in% rownames(spe_DLPFC_12),]

# drop rownames in loadings not in spe_DLPFC_30
loadings <- loadings[rownames(loadings) %in% rowData(spe_DLPFC_30)$gene_id,]

logcounts <- logcounts(spe_DLPFC_30)

proj <- project(loadings, as.matrix(logcounts))
proj <- t(proj)
colnames(proj) <- paste("NMF", 1:100, sep = "_")

# add to reducedDims
reducedDim(spe_DLPFC_30, "NMF_proj") <- proj

# save spe_DLPFC_30
save(spe_DLPFC_30, file = here("processed-data", "14_cross_region", "spe_DLPFC_30_NMF.Rdata"))

spe_DLPFC_30.temp <- spe_DLPFC_30

# add each proj column to colData(spe_DLPFC_30.temp)
for (i in 1:100){
    colData(spe_DLPFC_30.temp)[[paste0("NMF_",i)]] <- reducedDims(spe_DLPFC_30.temp)$NMF_proj[,i]
}

brains <- unique(spe_DLPFC_30.temp$subject)

for (i in 1:100){
    print(paste0("i=", i))

    pdf(file = here::here("plots", "14_cross_region", "SpotPlots_DLPFC_30", paste0("NMF_", i, ".pdf")),
        width = 21, height = 20)

    for (j in seq_along(brains)){
        spe_DLPFC_30b <- spe_DLPFC_30.temp[, which(spe_DLPFC_30.temp$subject == brains[j])]
        samples <- unique(spe_DLPFC_30b$sample_id)
        print(length(samples))

        if (length(samples) == 1){
            p1 <- vis_gene(spe =  spe_DLPFC_30b, sampleid = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 9, )
            grid.arrange(p1, nrow = 1)
        } else if (length(samples) == 2){
            p1 <- vis_gene(spe =  spe_DLPFC_30b, sampld = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p2 <- vis_gene(spe =  spe_DLPFC_30b, sampleid = samples[2], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            grid.arrange(p1, p2, nrow = 2)
        } else if (length(samples) == 3){
            p1 <- vis_gene(spe =  spe_DLPFC_30b, sampleid = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p2 <- vis_gene(spe =  spe_DLPFC_30b, sampleid = samples[2], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p3 <- vis_gene(spe =  spe_DLPFC_30b, sampleid = samples[3], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            grid.arrange(p1, p2, p3, nrow = 2)
        }
    }

    dev.off()
}

# find the number of zeroes in each column in colData(spe_DLPFC_30.temp)
zeroes <- sapply(colData(spe_DLPFC_30.temp), function(x) sum(x == 0))

# list out columns in colData(spe_DLPFC_30.temp) with more than 113917 zeroes
names(zeroes[zeroes > 113917])

