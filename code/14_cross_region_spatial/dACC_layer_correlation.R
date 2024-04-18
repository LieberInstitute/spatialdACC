setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(RcppML)
library(SpatialExperiment)
library(pheatmap)
library(reshape2)
library(here)

# load dACC NMF data
load(file = here("processed-data", "14_cross_region", "spe_dACC_NMF.Rdata"))

x <- reducedDim(spe_dACC, "NMF_proj")

# load dACC PRECAST data
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9.Rdata"))

spe_dACC$PRECAST_cluster <- spe$PRECAST_cluster

sum(is.na(spe_dACC$PRECAST_cluster))
# 0

spe_dACC$PRECAST_cluster <- unfactor(spe_dACC$PRECAST_cluster)
spe_dACC$PRECAST_cluster[spe_dACC$PRECAST_cluster == 3] <- "WM1"
spe_dACC$PRECAST_cluster[spe_dACC$PRECAST_cluster == 8] <- "WM2"
spe_dACC$PRECAST_cluster[spe_dACC$PRECAST_cluster == 7] <- "WM-CC"
spe_dACC$PRECAST_cluster[spe_dACC$PRECAST_cluster == 5] <- "L6b"
spe_dACC$PRECAST_cluster[spe_dACC$PRECAST_cluster == 6] <- "L6a"
spe_dACC$PRECAST_cluster[spe_dACC$PRECAST_cluster == 4] <- "L5"
spe_dACC$PRECAST_cluster[spe_dACC$PRECAST_cluster == 2] <- "L3"
spe_dACC$PRECAST_cluster[spe_dACC$PRECAST_cluster == 1] <- "L2"
spe_dACC$PRECAST_cluster[spe_dACC$PRECAST_cluster == 9] <- "L1"

####onehot encode layer
data<-as.data.frame(spe_dACC$PRECAST_cluster)
colnames(data)<-'layer'
onehot_layer <-  dcast(data = data, rownames(data) ~ layer, length)
rownames(onehot_layer)<-onehot_layer[,1]
onehot_layer[,1]<-as.numeric(onehot_layer[,1])
onehot_layer<-onehot_layer[order(onehot_layer[,1],decreasing=F),]
onehot_layer[,1]<-NULL

###correlate with nmf patterns
pdf(here("plots", "14_cross_region", "nmf_dACC_layer_correlation_heatmap.pdf"))
pheatmap(cor(x,onehot_layer), fontsize_row = 5)
dev.off()

# aggregate NMF patterns

# create dataframe
data <- data.frame(colData(spe_dACC), x)

# aggregate NMF patterns across layers
# grep "NMF" to get all NMF patterns
agg_data <- aggregate(data[,grep("NMF", colnames(data))],
                      by=list(data$PRECAST_cluster),
                      FUN=mean)

# move Group.1 to row names, then drop
rownames(agg_data) <- agg_data$Group.1
agg_data <- agg_data[,-1]

pdf(here("plots", "14_cross_region", "nmf_dACC_layer_correlation_aggregated_heatmap.pdf"))
p1 <- pheatmap(agg_data,
               color=colorRampPalette(c("blue","white","red"))(100),
               cluster_cols=T,
               cluster_rows=T,
               scale="column",
               fontsize_col = 5
)
dev.off()

