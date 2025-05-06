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

# save modeling results list
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
k <- k[k != "Sst Chodl"]

k1 <- k[c(1,3,4,9,12,15,17,18)]
k2 <- setdiff(k,k1)

plot_list <- list()

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

    p <- EnhancedVolcano(df_list[[i]],
                          lab = df_list[[i]]$gene_name,
                          pointSize = 1,
                          x = 'logFC',
                          y = 'FDR',
                          FCcutoff = 1.5,
                          pCutoff = 0.05,
                          ylab = "-log10 FDR",
                          legendLabels = c('Not sig.','Log (base 2) FC','FDR',
                                           'FDR & Log (base 2) FC'),
                         title = paste0(i, " vs. all others"),
                         subtitle = "",
                         caption = ""
    )

    plot_list[[i]] <- p
}

dev.off()


# supp figures
png(file = here::here("plots", "snRNA-seq","05_azimuth", "azimuth_DE_volcano_plots_layer.png"),
    width = 13, height = 18, unit="in", res=300)

wrap_plots(plot_list[["L2_3_IT"]],plot_list[["L5_ET"]],plot_list[["L5_IT"]],plot_list[["L5_6_NP"]],
           plot_list[["L6_CT"]],plot_list[["L6_IT"]],plot_list[["L6_IT_Car3"]],plot_list[["L6b"]],
           nrow=4, guides="collect") & theme(legend.position = 'bottom')

dev.off()


sort(k2)
png(file = here::here("plots", "snRNA-seq","05_azimuth", "azimuth_DE_volcano_plots_other.png"),
    width = 18, height = 18, unit="in", res=300)

wrap_plots(plot_list[["Astro"]],plot_list[["Endo"]],plot_list[["Lamp5"]],plot_list[["MicroPVM"]],
           plot_list[["Oligo"]],plot_list[["OPC"]],plot_list[["Pvalb"]],plot_list[["Sncg"]],
           plot_list[["Sst"]],plot_list[["Vip"]],plot_list[["VLMC"]],
           nrow=4, guides="collect") & theme(legend.position = 'bottom')

dev.off()

# DE between L5 ET and L5 IT
# want the positive x-axis to be L5 ET

# create volcano plot of pairwise comparison of NMF_38 and NMF_61
#volcano plots
thresh_fdr <- 0.05
thresh_logfc <- log2(1.5)

fdrs <- modeling_results[["pairwise"]][,paste0("fdr_", "L5_ET-L5_IT")]
logfc <- modeling_results[["pairwise"]][,paste0("logFC_", "L5_ET-L5_IT")]
#flip so L5 ET is positive
logfc <- -logfc

sig <- (fdrs < thresh_fdr) & (abs(logfc) > thresh_logfc)
print(table(sig))

df_list <- data.frame(
    gene_name = modeling_results[["pairwise"]]$gene,
    logFC = -logfc,
    FDR = fdrs,
    sig = sig
)

pdf(file = here::here("plots", "snRNA-seq", "05_azimuth", "volcano_plots_L5_IT_L5_ET.pdf"),
    width = 7, height = 6)

print(EnhancedVolcano(df_list,
                      lab = df_list$gene_name,
                      x = 'logFC',
                      y = 'FDR',
                      xlim = c(-3, 5),
                      ylim = c(0, -log10(10e-22)),
                      legendPosition = "bottom",
                      selectLab = c("VAT1L", "POU3F1", "SULF2", "HAPLN4", "LYPD1", "FEZF2", "GABRQ"),
                      FCcutoff = 1.5,
                      pCutoff = 0.05,
                      labSize = 7.0,
                      ylab = "-log10 FDR",
                      legendLabels = c('Not sig.','LogFC','FDR',
                                       'FDR & LogFC'),
                      title = "L5 ET vs. L5 IT",
                      subtitle = "",
                      caption = ""
)
)

dev.off()
