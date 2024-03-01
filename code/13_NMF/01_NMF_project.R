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

# get NMF results from single nucleus data
x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))

# load uncorrected and lognormalized spatial data
load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))

# load Single Nucleus object
# this file says "harmony", but that is only in the reduced dims, the counts/logcounts are not batch corrected
load(file = here("processed-data", "snRNA-seq", "03_batch_correction", "sce_harmony.Rdata"))

# extract patterns
patterns <- t(x$h)
colnames(patterns) <- paste("NMF", 1:120, sep = "_")

loadings <- x$w
rownames(loadings) <- rownames(sce)

# ====== project loadings to spatial data =======
# drop any gene_names in spe not in sce
# there were 168 genes in spe not in sce
spe <- spe[rowData(spe)$gene_name %in% rownames(sce),]

# drop 5314 rownames in loadings not in spe
loadings <- loadings[rownames(loadings) %in% rowData(spe)$gene_name,]

logcounts <- logcounts(spe)

proj <- project(logcounts, loadings)
proj <- t(proj)
colnames(proj) <- paste("NMF", 1:120, sep = "_")

# add to reducedDims
reducedDim(spe, "NMF_proj") <- proj

# save spe
save(spe, file = here("processed-data", "13_NMF", "spe_NMF.Rdata"))
