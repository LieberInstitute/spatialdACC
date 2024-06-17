setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(scuttle)
library(edgeR)

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))

# replace spaces with underscores
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L2/3 IT", "L2_3_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5 ET", "L5_ET")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5 IT", "L5_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5/6 NP", "L5_6_NP")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 CT", "L6_CT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 IT", "L6_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 IT Car3", "L6_IT_Car3")


sce<-SingleCellExperiment(assays=list(counts=counts(sce)),colData=colData(sce),rowData=rowData(sce))
sce_pseudo <- aggregateAcrossCells(
    sce,
    DataFrame(
        cluster = colData(sce)$cellType_azimuth
    ))

sce_pseudo$cluster <- factor(sce_pseudo$cluster)
sce_pseudo <- sce_pseudo[, sce_pseudo$ncells >= 50]
dim(sce_pseudo)
# [1] 34866    19

pdf(file = here::here("plots","17_LDSC", "snRNA-seq", "histogram_boxplot_azimuth.pdf"), width = 14, height = 14)
hist(sce_pseudo$ncells, breaks = 200)
boxplot(ncells ~ sce_pseudo$cluster, data = colData(sce_pseudo))
dev.off()

#find a good expression cutoff using edgeR::filterByExpr
rowData(sce_pseudo)$high_expr_group_cluster <- filterByExpr(sce_pseudo, group = sce_pseudo$cluster)

summary(rowData(sce_pseudo)$high_expr_group_cluster)
# Mode    FALSE    TRUE
# logical  6994   27872 

## Now filter
sce_pseudo <- sce_pseudo[rowData(sce_pseudo)$high_expr_group_cluster, ]
sce_pseudo <- sce_pseudo[!duplicated(rowData(sce_pseudo)$gene_name),]
sce_pseudo <- sce_pseudo[rowData(sce_pseudo)$gene_type=='protein_coding',]
dim(sce_pseudo)
# 16981    19

normd <- edgeR::cpm(edgeR::calcNormFactors(sce_pseudo))
rownames(normd)<-rowData(sce_pseudo)$gene_name
colnames(normd)<-sce_pseudo$cluster

write.table(normd, here::here("processed-data", "17_LDSC", "snRNA-seq_aggregated_cpm.tsv"), na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")

