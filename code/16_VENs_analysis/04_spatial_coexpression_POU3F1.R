setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(MERINGUE)
library(scuttle)

# load data
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata"))

# compute normcounts
spe <- logNormCounts(spe, transform = "none")

# 1. spatially unaware coexpression

# find index of VAT1L and POU3F1 in rowData(spe)$gene_name
VAT1L_index <- which(rowData(spe)$gene_name == "VAT1L")
POU3F1_index <- which(rowData(spe)$gene_name == "POU3F1")

dat <- data.frame(
    VAT1L = logcounts(spe)[VAT1L_index,],
    POU3F1 = logcounts(spe)[POU3F1_index,]
)

# compute general cross correlation
general_cor <- cor.test(dat$VAT1L, dat$POU3F1, method = "spearman")$estimate

# create scatterplot of gene expression VAT1L vs POU3F1
p <- ggplot(dat, aes(x = VAT1L, y = POU3F1)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    # print value of "general_cor" in the plot
    annotate("text", x = max(dat$VAT1L)-2, y = max(dat$POU3F1), label = paste("correlation:", round(general_cor, 2))) +
    labs(title = "VAT1L vs POU3F1 logcounts scale", x = "VAT1L expression", y = "POU3F1 expression") +
    theme_bw()

# 2. spatially aware coexpression
# separately by sample due to memory constraints

samples <- unique(colData(spe)$sample_id)

# save plots
pdf(here("plots", "16_VENs_analysis", "spatial_coexpression_POU3F1.pdf"))
print(p)

# create empty list to store results
results <- list()

for (sample in samples) {
    print(sample)
    spe_sample <- spe[, colData(spe)$sample_id == sample]
    pos <- data.frame(
        x = spe_sample$array_row,
        y = spe_sample$array_col
    )
    rownames(pos) <- colnames(counts(spe_sample))
    weight <- getSpatialNeighbors(pos,filterDist = 3)
    plotNetwork(pos, weight, cex = 0.5, main = sample)
    pval <- spatialCrossCorTest(normcounts(spe_sample)[VAT1L_index,], normcounts(spe_sample)[POU3F1_index,], weight,
                        ncores = 3, plot=T)
    print(pval)
    results[[sample]] <- pval
}

dev.off()

# save results as csv
write.csv(results, here("processed-data", "16_VENs_analysis", "spatial_coexpression_POU3F1.csv"))

