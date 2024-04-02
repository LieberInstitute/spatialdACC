setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

####code for correlating categorical technical variables with NMF patterns
library(here)
library(SpatialExperiment)
library(pheatmap)
library(reshape2)

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))
x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))

####onehot encode cell type
data<-as.data.frame(sce$cellType_azimuth)
colnames(data)<-'cellType_azimuth'
onehot_cellType_azimuth <-  dcast(data = data, rownames(data) ~ cellType_azimuth, length)
rownames(onehot_cellType_azimuth)<-onehot_cellType_azimuth[,1]
onehot_cellType_azimuth[,1]<-as.numeric(onehot_cellType_azimuth[,1])
onehot_cellType_azimuth<-onehot_cellType_azimuth[order(onehot_cellType_azimuth[,1],decreasing=F),]
onehot_cellType_azimuth[,1]<-NULL

###correlate with nmf patterns
pdf(here("plots", "snRNA-seq", "06_NMF", "nmf_cellType_azimuth_correlation_heatmap.pdf"))
pheatmap(cor(t(x@h),onehot_cellType_azimuth), fontsize_row = 5)
dev.off()

# aggregate NMF patterns

# create dataframe
data <- data.frame(colData(sce), t(x@h))

# aggregate NMF patterns across cell types
# grep "NMF" to get all NMF patterns
agg_data <- aggregate(data[,grep("nmf", colnames(data))],
                      by=list(data$cellType_azimuth),
                      FUN=mean)

# move Group.1 to row names, then drop
rownames(agg_data) <- agg_data$Group.1
agg_data <- agg_data[,-1]

pdf(here("plots", "snRNA-seq","06_NMF", "nmf_cellType_azimuth_correlation_aggregated_heatmap.pdf"))
p1 <- pheatmap(agg_data,
               color=colorRampPalette(c("blue","white","red"))(100),
               cluster_cols=T,
               cluster_rows=T,
               scale="column",
               fontsize_col = 5
)
dev.off()
