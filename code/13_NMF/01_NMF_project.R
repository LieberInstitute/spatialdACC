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
patterns <- t(x$h)
colnames(patterns) <- paste("NMF", 1:100, sep = "_")

loadings <- x$w
rownames(loadings) <- rownames(sce)

# ====== project loadings to spatial data =======
# drop any gene_names in spe not in sce
# there were 168 genes in spe not in sce
spe <- spe[rowData(spe)$gene_name %in% rownames(sce),]

# drop 5314 rownames in loadings not in spe
loadings <- loadings[rownames(loadings) %in% rowData(spe)$gene_name,]

logcounts <- logcounts(spe)

proj <- project(loadings, as.matrix(logcounts))
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

# find the number of zeroes in each column in colData(spe.temp)
zeroes <- sapply(colData(spe.temp)[, 1:100], function(x) sum(x == 0))

# list out columns in colData(spe.temp) with more than 77533 zeroes
names(zeroes[zeroes > 77533])
#  [4] "NMF_1"                "NMF_2"                "NMF_9"
# [7] "NMF_14"                "NMF_15"               "NMF_16"
# [10] "NMF_17"               "NMF_21"               "NMF_30"
# [13] "NMF_35"
