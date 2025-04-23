library(SingleCellExperiment)

load(file="normalized_dfs.Rdata")

# create single cell experiment
counts_mat <- cbind(t(Br6432_dACC_norm),t(Br6432_dlPFC_norm),
                    t(Br8325_dlPFC_norm),t(Br8325_left_dACC_norm),t(Br8325_right_dACC_norm))

row_mat <- data.frame(
    gene_name = rownames(counts_mat)
)

col_mat <- data.frame(
    sample = c(rep("Br6432_dACC",dim(Br6432_dACC_norm)[1]), rep("Br6432_dlPFC",dim(Br6432_dlPFC_norm)[1]), rep("Br8325_dlPFC",dim(Br8325_dlPFC_norm)[1]),
               rep("Br8325_left_dACC",dim(Br8325_left_dACC_norm)[1]), rep("Br8325_right_dACC",dim(Br8325_right_dACC_norm)[1])),
    region = c(rep("dACC",dim(Br6432_dACC_norm)[1]), rep("dlPFC",dim(Br6432_dlPFC_norm)[1]), rep("dlPFC",dim(Br8325_dlPFC_norm)[1]),
                     rep("dACC",dim(Br8325_left_dACC_norm)[1]), rep("dACC",dim(Br8325_right_dACC_norm)[1]))
)

sce <- SingleCellExperiment(list(counts=as.matrix(counts_mat)),
                            colData=col_mat,
                            rowData=row_mat)

# run k means clustering with 2 groups - hopefully VENs and not VENs
set.seed(1032)
res <- kmeans(t(counts(sce)),
              centers = 3)
colData(sce)$cluster <- res$cluster

table(colData(sce)$cluster, colData(sce)$region)
table(colData(sce)$cluster, colData(sce)$sample)

pca <- prcomp(t(assays(sce)$counts))
metadata(sce) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
pca_pseudo <- pca$x[, seq_len(4)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(sce) <- list(PCA = pca_pseudo)

# plot counts of clusters from dACC or dlPFC
p1 <- plotPCA(
    sce,
    colour_by = "region",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(sce)$PCA_var_explained
)

p2 <- plotPCA(
    sce,
    colour_by = "sample",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(sce)$PCA_var_explained
)

p3 <- plotPCA(
    sce,
    colour_by = "cluster",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(sce)$PCA_var_explained
)

wrap_plots(p1,p2,p3)

# look at gene loadings to see what is driving this
pcar <- pca$rotation
barplot(pcar[, c('PC1', 'PC2')], beside=TRUE, col=seq_len(ncol(pcar)) + 1,
        legend=rownames(pcar), args.legend=list(x='top'))
abline(h=0.3, col='red', lty=2)

