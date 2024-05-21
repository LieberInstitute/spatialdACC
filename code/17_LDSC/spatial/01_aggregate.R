setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(scuttle)
library(edgeR)

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


spe<-SingleCellExperiment(assays=list(counts=counts(spe)),colData=colData(spe),rowData=rowData(spe))
spe_pseudo <- aggregateAcrossCells(
    spe,
    DataFrame(
        cluster = colData(spe)$PRECAST_cluster
    ))

spe_pseudo$cluster <- factor(spe_pseudo$cluster)
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 50]
dim(spe_pseudo)
# [1] 29720     9

pdf(file = here::here("plots","17_LDSC", "spatial", "histogram_boxplot_precast9.pdf"), width = 14, height = 14)
hist(spe_pseudo$ncells, breaks = 200)
boxplot(ncells ~ spe_pseudo$cluster, data = colData(spe_pseudo))
dev.off()

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$cluster)

summary(rowData(spe_pseudo)$high_expr_group_cluster)
# Mode    FALSE    TRUE
# logical 10560   19160

## Now filter
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_cluster, ]
spe_pseudo <- spe_pseudo[!duplicated(rowData(spe_pseudo)$gene_name),]
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$gene_type=='protein_coding',]
dim(spe_pseudo)
# 15293     9

normd <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo))
rownames(normd)<-rowData(spe_pseudo)$gene_name
colnames(normd)<-spe_pseudo$cluster

write.table(normd, here::here("processed-data", "17_LDSC", "visium_aggregated_cpm.tsv"), na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")

