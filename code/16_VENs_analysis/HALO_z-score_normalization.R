library(ggplot2)
library(tidyverse)
library(patchwork)
library(preprocessCore)
library(SingleCellExperiment)
library(factoextra)
library(cluster)
library(mbkmeans)
library(scater)

setwd("../../../Downloads/Data")

load(file="nucleus_normalized_dfs.Rdata")

# calculate z-scores for each gene within each sample

Br6432_dACC_norm <- as.data.frame(scale(Br6432_dACC))
Br6432_dlPFC_norm <- as.data.frame(scale(Br6432_dlPFC))
Br8325_dlPFC_norm <- as.data.frame(scale(Br8325_dlPFC))
Br8325_left_dACC_norm <- as.data.frame(scale(Br8325_left_dACC))
Br8325_right_dACC_norm <- as.data.frame(scale(Br8325_right_dACC))

# vis boxplots of expression after z-score normalization
datasets <- list(
    Br6432_dACC = Br6432_dACC_norm,
    Br6432_dlPFC = Br6432_dlPFC_norm,
    Br8325_dlPFC = Br8325_dlPFC_norm,
    Br8325_left_dACC = Br8325_left_dACC_norm,
    Br8325_right_dACC = Br8325_right_dACC_norm
)

long_df <- bind_rows(
    lapply(names(datasets), function(name) {
        datasets[[name]] %>%
            mutate(dataset = name) %>%
            pivot_longer(
                cols = -dataset,
                names_to = "gene",
                values_to = "expression"
            )
    })
)

long_df$expression <- long_df$expression +1

ggplot(long_df, aes(x = gene, y = expression, color=dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression (log10(x+1) scale)", color = "Dataset") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(trans='log10')

long_df$expression <- long_df$expression +1

p1 <- ggplot(long_df[which(long_df$gene=="GABRQ"),], aes(x = gene, y = expression, color = dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression", color = "Dataset") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(trans="log2")
p2 <- ggplot(long_df[which(long_df$gene=="PCP4"),], aes(x = gene, y = expression, color = dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression", color = "Dataset") +
    scale_fill_brewer(palette = "Set2")+
    scale_y_continuous(trans="log2")
p3 <- ggplot(long_df[which(long_df$gene=="POU3F1"),], aes(x = gene, y = expression, color = dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression", color = "Dataset") +
    scale_fill_brewer(palette = "Set2")+
    scale_y_continuous(trans="log2")
p4 <- ggplot(long_df[which(long_df$gene=="SULF2"),], aes(x = gene, y = expression, color = dataset)) +
    geom_boxplot(position = position_dodge(0.8)) +
    theme_minimal() +
    labs(x = "Gene", y = "Expression", color = "Dataset") +
    scale_fill_brewer(palette = "Set2")+
    scale_y_continuous(trans="log2")

wrap_plots(p1,p2,p3,p4, nrow=2, guides="collect")

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

pca <- prcomp(t(assays(sce)$counts))
metadata(sce) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
pca_pseudo <- pca$x[, seq_len(4)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(sce) <- list(PCA = pca_pseudo)

# look at gene loadings to see what is driving this
pcar <- pca$rotation
barplot(pcar[, c('PC1', 'PC2')], beside=TRUE, col=seq_len(ncol(pcar)) + 1,
        legend=rownames(pcar), args.legend=list(x='top'))
abline(h=0.3, col='red', lty=2)

# look at elbow plot
ks <- seq(2, 15)
res <- lapply(ks, function(k) {
    mbkmeans(sce, clusters = k,
             reduceMethod = NA,
             calc_wcss = TRUE, num_init=10)
})
wcss <- sapply(res, function(x) sum(x$WCSS_per_cluster))
plot(ks, wcss, type = "b")

# run k means clustering with 2 groups - hopefully VENs and not VENs
set.seed(1032)
res <- kmeans(t(counts(sce)),
              centers = 5)
colData(sce)$cluster <- res$cluster
colData(sce)$cluster <- as.factor(colData(sce)$cluster)

table(colData(sce)$cluster, colData(sce)$region)
table(colData(sce)$cluster, colData(sce)$sample)

plotExpression(sce, features=rowData(sce)$gene_name, x="cluster",colour_by = "cluster", exprs_values="counts")


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

