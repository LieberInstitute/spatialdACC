setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(scuttle)
library(edgeR)

load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata"))

spe<-SingleCellExperiment(assays=list(counts=counts(spe)),colData=colData(spe),rowData=rowData(spe))
spe_pseudo <- aggregateAcrossCells(
    spe,
    DataFrame(
        cluster = colData(spe)$layer
    ))

spe_pseudo$cluster <- factor(spe_pseudo$cluster)
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 50]
dim(spe_pseudo)
# [1] 29720     7

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$cluster)

summary(rowData(spe_pseudo)$high_expr_group_cluster)

## Now filter
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_cluster, ]
spe_pseudo <- spe_pseudo[!duplicated(rowData(spe_pseudo)$gene_name),]
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$gene_type=='protein_coding',]
dim(spe_pseudo)
# 15814     7

normd <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo))
rownames(normd)<-rowData(spe_pseudo)$gene_name
colnames(normd)<-spe_pseudo$cluster

write.table(normd, here::here("processed-data", "17_LDSC", "visium_aggregated_cpm.tsv"), na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")

