setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SingleCellExperiment)
library(spatialLIBD)
library(here)
library(scran)
library('EnhancedVolcano')

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))
sce <- logNormCounts(sce)

# replace spaces with underscores
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L2/3 IT", "L2_3_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5 ET", "L5_ET")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5 IT", "L5_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5/6 NP", "L5_6_NP")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 CT", "L6_CT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 IT", "L6_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 IT Car3", "L6_IT_Car3")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "Sst Chodl", "Sst_Chodl")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "Micro-PVM", "Micro_PVM")

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

sce_pseudo <-
    registration_pseudobulk(sce,
                            var_registration = "cellType_azimuth",
                            var_sample_id = "brain"
    )

sig_genes <- sig_genes_extract(
    n = 20,
    modeling_results = modeling_results,
    model_type = "enrichment",
    sce_layer = sce_pseudo
)

save(
    sig_genes,
    file = here("processed-data", "snRNA-seq", "05_azimuth", "azimuth_DE_sig_genes_top20.Rdata")
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

