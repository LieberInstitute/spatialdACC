setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

####code for correlating PRECAST cluster with NMF patterns
library(here)
library(SpatialExperiment)
library(pheatmap)
library(reshape2)

load(file = here("processed-data", "13_NMF", "spe_NMF.Rdata"))

# proj is the columns in colData(spe) that start with "nmf"
proj <- as.matrix(colData(spe)[,grep("nmf", colnames(colData(spe)))])

# load precast cluster data
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata"))

spe$PRECAST_cluster <- spe$layer

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

# aggregate NMF patterns

# create dataframe
data <- data.frame(colData(spe), proj_no_zero)

# aggregate NMF patterns across clusters. # grep "NMF" to get all NMF patterns
agg_data <- aggregate(data[,grep("nmf", colnames(data))],
                      by=list(data$PRECAST_cluster),
                      FUN=mean)

# move Group.1 to row names, then drop
rownames(agg_data) <- agg_data$Group.1
agg_data <- agg_data[,-1]

pdf(here("plots", "13_NMF", "nmf_precast_correlation_aggregated_heatmap.pdf"))
p1 <- pheatmap(agg_data,
               color=colorRampPalette(c("blue","white","red"))(100),
               cluster_cols=T,
               cluster_rows=T,
               scale="column",
               fontsize_col = 5
)
dev.off()

####onehot encode sample_id
data<-as.data.frame(spe$sample_id)
colnames(data)<-'sample_id'
onehot_sample_id <-  dcast(data = data, rownames(data) ~ sample_id, length)
rownames(onehot_sample_id)<-onehot_sample_id[,1]
onehot_sample_id[,1]<-as.numeric(onehot_sample_id[,1])
onehot_sample_id<-onehot_sample_id[order(onehot_sample_id[,1],decreasing=F),]
onehot_sample_id[,1]<-NULL

####onehot encode slide
data<-as.data.frame(spe$slide)
colnames(data)<-'slide'
onehot_slide <-  dcast(data = data, rownames(data) ~ slide, length)
rownames(onehot_slide)<-onehot_slide[,1]
onehot_slide[,1]<-as.numeric(onehot_slide[,1])
onehot_slide<-onehot_slide[order(onehot_slide[,1],decreasing=F),]
onehot_slide[,1]<-NULL

####onehot encode brnum
data<-as.data.frame(spe$brnum)
colnames(data)<-'brnum'
onehot_brnum <-  dcast(data = data, rownames(data) ~ brnum, length)
rownames(onehot_brnum)<-onehot_brnum[,1]
onehot_brnum[,1]<-as.numeric(onehot_brnum[,1])
onehot_brnum<-onehot_brnum[order(onehot_brnum[,1],decreasing=F),]
onehot_brnum[,1]<-NULL

###correlate with nmf patterns
pdf(here("plots", "13_NMF", "nmf_sample_correlation_heatmap.pdf"))
pheatmap(cor(proj_no_zero,as.matrix(onehot_sample_id)), fontsize_row = 5)
pheatmap(cor(proj_no_zero,as.matrix(onehot_slide)), fontsize_row = 5)
pheatmap(cor(proj_no_zero,as.matrix(onehot_brnum)), fontsize_row = 5)
dev.off()

