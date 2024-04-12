setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

####code for correlating categorical technical variables with NMF patterns
library(here)
library(SpatialExperiment)
library(pheatmap)
library(reshape2)

# get NMF results from DLPFC data
x <- readRDS(file = here("processed-data", "14_cross_region", "DLPFC_nmf_results.RDS"))
# load DLPFC spatial object
spe_DLPFC <- spatialLIBD::fetch_data(type = "spe")

sum(is.na(spe_DLPFC$layer_guess_reordered))
#[1] 352

# remove NAs
spe_DLPFC <- spe_DLPFC[,!is.na(spe_DLPFC$layer_guess_reordered)]

# remove NAs from NMF results
patterns <- x@h
patterns <- patterns[,!is.na(spe_DLPFC$layer_guess_reordered)]

####onehot encode layer
data<-as.data.frame(spe_DLPFC$layer_guess_reordered)
colnames(data)<-'layer'
onehot_layer <-  dcast(data = data, rownames(data) ~ layer, length)
rownames(onehot_layer)<-onehot_layer[,1]
onehot_layer[,1]<-as.numeric(onehot_layer[,1])
onehot_layer<-onehot_layer[order(onehot_layer[,1],decreasing=F),]
onehot_layer[,1]<-NULL

###correlate with nmf patterns
pdf(here("plots", "14_cross_region", "nmf_DLPFC_layer_correlation_heatmap.pdf"))
pheatmap(cor(t(patterns),onehot_layer), fontsize_row = 5)
dev.off()

# aggregate NMF patterns

# create dataframe
data <- data.frame(colData(spe_DLPFC), t(patterns))

# aggregate NMF patterns across layers
# grep "NMF" to get all NMF patterns
agg_data <- aggregate(data[,grep("nmf", colnames(data))],
                      by=list(data$layer_guess_reordered),
                      FUN=mean)

# move Group.1 to row names, then drop
rownames(agg_data) <- agg_data$Group.1
agg_data <- agg_data[,-1]

pdf(here("plots", "14_cross_region", "nmf_DLPFC_layer_correlation_aggregated_heatmap.pdf"))
p1 <- pheatmap(agg_data,
               color=colorRampPalette(c("blue","white","red"))(100),
               cluster_cols=T,
               cluster_rows=T,
               scale="column",
               fontsize_col = 5
)
dev.off()
