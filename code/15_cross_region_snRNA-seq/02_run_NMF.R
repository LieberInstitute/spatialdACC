setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(RcppML)
library(pheatmap)
library(SingleCellExperiment)
library(here)
library(scuttle)
library(Matrix)
library(spatialLIBD)
library(reshape2)

# Load data
sce_path_zip <- fetch_data("spatialDLPFC_snRNAseq")
sce_path <- unzip(sce_path_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
    file.path(tempdir(), "sce_DLPFC_annotated")
)

assay(sce, "logcounts") <- as(assay(sce, "logcounts"), "dgCMatrix")

options(RcppML.threads=4)
x <- RcppML::nmf(assay(sce,'logcounts'),
                 k=75,
                 tol = 1e-06,
                 maxit = 1000,
                 verbose = T,
                 L1 = 0.1,
                 seed = 135,
                 mask_zeros = FALSE,
                 diag = TRUE,
                 nonneg = TRUE
)

# Save NMF results
saveRDS(x, file = here("processed-data", "15_cross_region_snRNA-seq", "DLPFC_nmf_results.RDS"))

####onehot encode layer
data<-as.data.frame(sce$layer_annotation)
colnames(data)<-'layer'
onehot_layer <-  dcast(data = data, rownames(data) ~ layer, length)
rownames(onehot_layer)<-onehot_layer[,1]
onehot_layer[,1]<-as.numeric(onehot_layer[,1])
onehot_layer<-onehot_layer[order(onehot_layer[,1],decreasing=F),]
onehot_layer[,1]<-NULL

###correlate with nmf patterns
pdf(here("plots", "15_cross_region_snRNA-seq", "nmf_DLPFC_layer_correlation_heatmap.pdf"))
pheatmap(cor(t(x@h),onehot_layer), fontsize_row = 5)
dev.off()

# aggregate NMF patterns

# create dataframe
data <- data.frame(colData(sce), t(x@h))

# aggregate NMF patterns across layers
# grep "NMF" to get all NMF patterns
agg_data <- aggregate(data[,grep("nmf", colnames(data))],
                      by=list(data$layer_annotation),
                      FUN=mean)

# move Group.1 to row names, then drop
rownames(agg_data) <- agg_data$Group.1
agg_data <- agg_data[,-1]

pdf(here("plots", "15_cross_region_snRNA-seq", "nmf_DLPFC_layer_correlation_aggregated_heatmap.pdf"))
p1 <- pheatmap(agg_data,
               color=colorRampPalette(c("blue","white","red"))(100),
               cluster_cols=T,
               cluster_rows=T,
               scale="column",
               fontsize_col = 5
)
dev.off()


pdf(here("plots", "15_cross_region_snRNA-seq", "nmf_DLPFC_layer_correlation_aggregated_unscaled_heatmap.pdf"))
p1 <- pheatmap(agg_data,
               color=colorRampPalette(c("blue","white","red"))(100),
               cluster_cols=T,
               cluster_rows=T,
               scale="none",
               fontsize_col = 5
)
dev.off()

pdf(here("plots", "15_cross_region_snRNA-seq", "nmf_DLPFC_layer_correlation_aggregated_rowscale_heatmap.pdf"))
p1 <- pheatmap(agg_data,
               color=colorRampPalette(c("blue","white","red"))(100),
               cluster_cols=T,
               cluster_rows=T,
               scale="row",
               fontsize_col = 5
)
dev.off()

nmf_df <- as.data.frame(t(x$h)) # the patterns have the correct dims, but in the comment it says loadings

# normalize loadings
col_sums <- colSums(nmf_df) # these are all = 1?

# Normalize each column by its sum
normalized_df <- sweep(nmf_df, 2, col_sums, FUN="/")
normalized_df$group <- sce$layer_annotation
normalized_df <- normalized_df[!is.na(normalized_df) != 0,] # this makes the dims get large

summarized <- summarizeAssayByGroup(t(normalized_df),
                                    ids=DataFrame(group="group"),
                                    #subset.row=subset_patterns,
                                    statistics=c("sum", "prop.detected"), threshold=0) # this does not work

sum <- assay(summarized, "sum")
num <- assay(summarized, "prop.detected")
group.names <- summarized$group

evals_long <- data.frame(
    Feature=rep(colnames(normalized_df), ncol(num)),
    Group=rep(group.names, each=nrow(num)),
    NumDetected=as.numeric(num),
    Sum=as.numeric(sum)
)


