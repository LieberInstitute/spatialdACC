setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(MERINGUE)
library(scuttle)

# load data
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9.Rdata"))

# compute normcounts
spe <- logNormCounts(spe, transform = "none")

# 1. spatially unaware coexpression

# find index of VAT1L and DRD5 in rowData(spe)$gene_name
VAT1L_index <- which(rowData(spe)$gene_name == "VAT1L")
DRD5_index <- which(rowData(spe)$gene_name == "DRD5")

dat <- data.frame(
    VAT1L = normcounts(spe)[VAT1L_index,],
    DRD5 = normcounts(spe)[DRD5_index,]
)

# create scatterplot of gene expression VAT1L vs DRD5
p <- ggplot(dat, aes(x = VAT1L, y = DRD5)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = "VAT1L vs DRD5 normcounts scale", x = "VAT1L", y = "DRD5")

# compute general cross correlation
cor.test(dat$VAT1L, dat$DRD5, method = "spearman")

# 2. spatially aware coexpression

pos <- data.frame(
    x = spe$array_row,
    y = spe$array_col
)

weight <- getSpatialNeighbors(pos, filterDist = 1)
spatialCrossCor(dat$VAT1L, dat$DRD5, weight)

# save plots
pdf(here("plots", "16_VENs_analysis", "spatial_coexpression.pdf"))
print(p)
plotNetwork(pos, weight)
spatialCrossCorTest(dat$VAT1L, dat$DRD5, weight,
                    plot=TRUE)
dev.off()

