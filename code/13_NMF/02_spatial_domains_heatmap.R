setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

####code for correlating PRECAST cluster with NMF patterns
library(here)
library(SpatialExperiment)
library(pheatmap)
library(reshape2)

load(file = here("processed-data", "13_NMF", "spe_NMF.Rdata"))
proj <- reducedDim(spe, "NMF_proj")

# load precast cluster data
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9.Rdata"))

####onehot encode precast cluster
data<-as.data.frame(spe$PRECAST_cluster)
colnames(data)<-'precast'
onehot_precast <-  dcast(data = data, rownames(data) ~ precast, length)
rownames(onehot_precast)<-onehot_precast[,1]
onehot_precast[,1]<-as.numeric(onehot_precast[,1])
onehot_precast<-onehot_precast[order(onehot_precast[,1],decreasing=F),]
onehot_precast[,1]<-NULL

# fix error In cor(proj, onehot_precast) : the standard deviation is zero
# check if any of the columns in proj are all 0
sum(colSums(proj) == 0)

# remove columns in proj that are all 0
proj_no_zero <- proj[, colSums(proj) != 0]

###correlate with nmf patterns
pdf(here("plots", "13_NMF", "nmf_precast_correlation_heatmap.pdf"))
pheatmap(cor(proj_no_zero,as.matrix(onehot_precast)), fontsize_row = 5)
dev.off()

