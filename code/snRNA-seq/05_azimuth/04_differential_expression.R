setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SingleCellExperiment)
library(spatialLIBD)
library(here)
library(scran)
library('EnhancedVolcano')
library("dplyr")
library("scater")
library("patchwork")

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))
sum_by_sample <- setNames(aggregate(sum ~ Sample, colData(sce), sum), c("Sample", "sum_sample"))
detected_by_sample <- setNames(aggregate(detected ~ Sample, colData(sce), sum), c("Sample", "detected_sample"))

sce_pseudo <-
    registration_pseudobulk(sce,
                            var_registration = "cellType_azimuth",
                            var_sample_id = "brain"
    )

save(sce_pseudo, file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth_pseudo.Rdata"))

pca <- prcomp(t(assays(sce_pseudo)$logcounts))
metadata(sce_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
pca_pseudo <- pca$x[, seq_len(20)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(sce_pseudo) <- list(PCA = pca_pseudo)

## Plot PCs
col_data_df <- as.data.frame(colData(sce_pseudo))
col_data_df <- left_join(col_data_df, detected_by_sample, by = "Sample")
col_data_df <- left_join(col_data_df, sum_by_sample, by = "Sample")

colData(sce_pseudo) <- DataFrame(col_data_df)

# make supp figure
p1 <- plotPCA(
    sce_pseudo,
    colour_by = "cellType_azimuth",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(sce_pseudo)$PCA_var_explained
)

p2 <- plotPCA(
    sce_pseudo,
    colour_by = "brain",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(sce_pseudo)$PCA_var_explained
)

p3 <- plotPCA(
    sce_pseudo,
    colour_by = "sum_sample",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(sce_pseudo)$PCA_var_explained
)

vars <- getVarianceExplained(sce_pseudo,
                             variables = c("cellType_azimuth","brain", "sum_sample", "detected_sample")
)

p4 <- plotExplanatoryVariables(vars)

png(file = here("plots", "snRNA-seq",
                "05_azimuth",
                "pseudobulk_PC_azimuth.png"),
    width = 10, height = 10, unit="in", res=300)

wrap_plots(p1,p2,p3,p4,nrow=2)
dev.off()


modeling_results <- registration_wrapper(
    sce,
    var_registration = "cellType_azimuth",
    var_sample_id = "brain",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

###save modeling results list
save(
    modeling_results,
    file = here("processed-data", "snRNA-seq", "05_azimuth", "azimuth_DE_results.Rdata")
)

sig_genes <- sig_genes_extract(
    n = 30,
    modeling_results = modeling_results,
    model_type = "enrichment",
    sce_layer = sce_pseudo
)

save(
    sig_genes,
    file = here("processed-data", "snRNA-seq", "05_azimuth", "azimuth_DE_sig_genes_top30.Rdata")
)

write.csv(sig_genes, file = here("processed-data", "snRNA-seq", "05_azimuth", "azimuth_DE_sig_genes_top30.csv"),
          row.names = F
)


#volcano plots
thresh_fdr <- 0.05
thresh_logfc <- log2(1.5)
fdrs_gene_ids <- rowData(sce_pseudo)$gene_id
fdrs_gene_names <- rowData(sce_pseudo)$gene_name

df_list <- list()
k <- unique(sce$cellType_azimuth)
k <- k[k != "Sst_Chodl"]

pdf(file = here::here("plots", "snRNA-seq","05_azimuth", "azimuth_DE_volcano_plots.pdf"),
    width = 8.5, height = 8)

for (i in k) {

    print(i)
    fdrs <- modeling_results[["enrichment"]][,paste0("fdr_", i)]
    logfc <- modeling_results[["enrichment"]][,paste0("logFC_", i)]

    # Identify significant genes (low FDR and high logFC)
    sig <- (fdrs < thresh_fdr) & (abs(logfc) > thresh_logfc)

    # Number of significant genes
    print(paste("Cluster", i))
    print(table(sig))

    df_list[[i]] <- data.frame(
        gene_name = modeling_results[["enrichment"]]$gene,
        logFC = logfc,
        FDR = fdrs,
        sig = sig
    )

    print(EnhancedVolcano(df_list[[i]],
                          lab = df_list[[i]]$gene_name,
                          x = 'logFC',
                          y = 'FDR',
                          FCcutoff = 1.5,
                          pCutoff = 0.05,
                          ylab = "-log10 FDR",
                          legendLabels = c('Not sig.','Log (base 2) FC','FDR',
                                           'FDR & Log (base 2) FC'),
                          title = "Cell Type (Azimuth) dACC",
                          subtitle = paste0("Cluster ", i, " vs. all others")
    ))
}

dev.off()

