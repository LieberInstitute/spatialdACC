setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(RcppML)
library(pheatmap)
library(SingleCellExperiment)
library(here)
library(scuttle)
library(Matrix)
library(reshape2)
library(spatialLIBD)

# Load data
spe <- spatialLIBD::fetch_data(type = "spe")

# want to use logcounts for NMF
spe <- logNormCounts(spe)

options(RcppML.threads=4)
x <- RcppML::nmf(assay(spe,'logcounts'),
		 k=100,
                 tol = 1e-06,
                 maxit = 1000,
                 verbose = T,
                 L1 = 0.1,
                 seed = 1135,
                 mask_zeros = FALSE,
                 diag = TRUE,
                 nonneg = TRUE
)

# Save NMF results
saveRDS(x, file = here("processed-data", "14_cross_region", "DLPFC_nmf_results.RDS"))

####onehot encode layer
data<-as.data.frame(spe$layer_guess_reordered)
colnames(data)<-'layer'
onehot_layer <-  dcast(data = data, rownames(data) ~ layer, length)
rownames(onehot_layer)<-onehot_layer[,1]
onehot_layer[,1]<-as.numeric(onehot_layer[,1])
onehot_layer<-onehot_layer[order(onehot_layer[,1],decreasing=F),]
onehot_layer[,1]<-NULL

###correlate with nmf patterns
pdf(here("plots", "14_cross_region", "nmf_DLPFC_12_layer_correlation_heatmap.pdf"))
pheatmap(cor(t(x@h),onehot_layer), fontsize_row = 5)
dev.off()

# aggregate NMF patterns

# create dataframe
data <- data.frame(colData(spe), t(x@h))

# aggregate NMF patterns across layers
# grep "NMF" to get all NMF patterns
agg_data <- aggregate(data[,grep("nmf", colnames(data))],
                      by=list(data$layer_guess_reordered),
                      FUN=mean)

# move Group.1 to row names, then drop
rownames(agg_data) <- agg_data$Group.1
agg_data <- agg_data[,-1]

pdf(here("plots", "14_cross_region", "nmf_DLPFC_12_layer_correlation_aggregated_heatmap.pdf"))
p1 <- pheatmap(agg_data,
               color=colorRampPalette(c("blue","white","red"))(100),
               cluster_cols=T,
               cluster_rows=T,
               scale="column",
               fontsize_col = 5
)
dev.off()
