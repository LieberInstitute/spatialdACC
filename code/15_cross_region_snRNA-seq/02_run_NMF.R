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

####onehot encode celltype
data<-as.data.frame(sce$cellType_hc)
colnames(data)<-'celltype'
onehot_celltype <-  dcast(data = data, rownames(data) ~ celltype, length)
rownames(onehot_celltype)<-onehot_celltype[,1]
onehot_celltype[,1]<-as.numeric(onehot_celltype[,1])
onehot_celltype<-onehot_celltype[order(onehot_celltype[,1],decreasing=F),]
onehot_celltype[,1]<-NULL

###correlate with nmf patterns
pdf(here("plots", "15_cross_region_snRNA-seq", "nmf_DLPFC_celltype_correlation_heatmap.pdf"))
pheatmap(cor(t(x@h),onehot_celltype), fontsize_row = 5)
dev.off()

# aggregate NMF patterns

# create dataframe
data <- data.frame(colData(sce), t(x@h))

# aggregate NMF patterns across celltypes
# grep "NMF" to get all NMF patterns
agg_data <- aggregate(data[,grep("nmf", colnames(data))],
                      by=list(data$cellType_hc),
                      FUN=mean)

# move Group.1 to row names, then drop
rownames(agg_data) <- agg_data$Group.1
agg_data <- agg_data[,-1]

pdf(here("plots", "15_cross_region_snRNA-seq", "nmf_DLPFC_celltype_correlation_aggregated_heatmap.pdf"))
p1 <- pheatmap(agg_data,
               color=colorRampPalette(c("blue","white","red"))(100),
               cluster_cols=T,
               cluster_rows=T,
               scale="column",
               fontsize_col = 5
)
dev.off()

