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
x <- readRDS(file = here("processed-data", "15_cross_region_snRNA-seq", "DLPFC_nmf_results.RDS"))

# load old DLPFC spatial object
spe_DLPFC_12 <- spatialLIBD::fetch_data(type = "spe")

# load DLPFC sce object
sce_path_zip <- fetch_data("spatialDLPFC_snRNAseq")
sce_path <- unzip(sce_path_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
    file.path(tempdir(), "sce_DLPFC_annotated")
)

# extract patterns
patterns <- t(x$h)
colnames(patterns) <- paste("NMF", 1:75, sep = "_")

loadings <- x$w
rownames(loadings) <- rowData(sce)$gene_id

# ====== project loadings to spatial data =======
# drop any gene_names in spe_DLPFC_12 not in speLPFC
spe_DLPFC_12 <- spe_DLPFC_12[rowData(spe_DLPFC_12)$gene_id %in% rowData(sce)$gene_id,]

# drop rownames in loadings not in spe_DLPFC_12
loadings <- loadings[rownames(loadings) %in% rowData(spe_DLPFC_12)$gene_id,]

logcounts <- logcounts(spe_DLPFC_12)

dim(loadings)
dim(logcounts)

proj <- project(loadings, as.matrix(logcounts))
proj <- t(proj)
colnames(proj) <- paste("NMF", 1:75, sep = "_")

# add to reducedDims
reducedDim(spe_DLPFC_12, "NMF_proj") <- proj

# save spe_DLPFC_12
save(spe_DLPFC_12, file = here("processed-data", "15_cross_region_snRNA-seq", "spe_DLPFC_12_NMF.Rdata"))

spe_DLPFC_12.temp <- spe_DLPFC_12

# add each proj column to colData(spe_DLPFC_12.temp)
for (i in 1:75){
    colData(spe_DLPFC_12.temp)[[paste0("NMF_",i)]] <- reducedDims(spe_DLPFC_12.temp)$NMF_proj[,i]
}

brains <- unique(spe_DLPFC_12.temp$subject)

for (i in 1:75){
    print(paste0("i=", i))

    pdf(file = here::here("plots", "15_cross_region_snRNA-seq", "SpotPlots_DLPFC_12", paste0("NMF_", i, ".pdf")),
        width = 21, height = 20)

    for (j in seq_along(brains)){
        spe_DLPFC_12b <- spe_DLPFC_12.temp[, which(spe_DLPFC_12.temp$subject == brains[j])]
        samples <- unique(spe_DLPFC_12b$sample_id)
        print(length(samples))

        if (length(samples) == 1){
            p1 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 9, )
            grid.arrange(p1, nrow = 1)
        } else if (length(samples) == 2){
            p1 <- vis_gene(spe =  spe_DLPFC_12b, sampld = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p2 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[2], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            grid.arrange(p1, p2, nrow = 2)
        } else if (length(samples) == 3){
            p1 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p2 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[2], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p3 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[3], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            grid.arrange(p1, p2, p3, nrow = 2)
        }
	else if (length(samples) == 4){
            p1 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p2 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[2], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p3 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[3], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p4 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[4], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            grid.arrange(p1, p2, p3, p4, nrow = 2)
        }
    }

    dev.off()
}

# find the number of zeroes in each column in colData(spe_DLPFC_12.temp)
zeroes <- sapply(colData(spe_DLPFC_12.temp)[, 1:75], function(x) sum(x == 0))

# list out columns in colData(spe_DLPFC_12.temp) with more than 77533 zeroes
names(zeroes[zeroes > 77533])

