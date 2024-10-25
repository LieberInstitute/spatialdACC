setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("scater")
    library("spatialLIBD")
    library("dplyr")
    library('EnhancedVolcano')
})

# read_barcoded_csv <- function(x) {
#     df <- read.csv(x)
#     colnames(df) <- tolower(colnames(df))
#
#     if (colnames(df)[2] == "cluster") {
#         colnames(df)[2] <-
#             gsub("gene_expression_", "", basename(dirname(x)))
#     }
#     return(df)
# }
#
# nnSVG_PRECAST_import <- function(spe, cluster_dir = file.path(tempdir(), "exported_clusters")) {
#     clustering_files <-
#         list.files(
#             here::here("processed-data", "08_clustering", "PRECAST", nnSVG_precast_name),
#             pattern = "clusters.csv",
#             all.files = TRUE,
#             full.names = TRUE,
#             recursive = TRUE
#         )
#     clusters_list <- lapply(clustering_files, read_barcoded_csv)
#     clusters <- Reduce(function(...) merge(..., by = "key", all = TRUE), clusters_list)
#     cluster_cols <- which(colnames(clusters) != "key")
#     colnames(clusters)[cluster_cols] <- paste0("", colnames(clusters)[cluster_cols])
#
#     colData(spe)[,nnSVG_precast_name] <- clusters[,2]
#     return(spe)
# }
#
# load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))
#
# k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#
# #import precast nnSVG clusters
# #i re wrote the key when using PRECAST, just add the cluster columns manually
# nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", k)
# spe <- nnSVG_PRECAST_import(
#     spe,
#     cluster_dir = here::here("processed-data", "08_clustering", "PRECAST", nnSVG_precast_name)
# )


# load spe for k=9 without WM-CC
load(
    file = here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata")
)

nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", 9)
colData(spe)[,nnSVG_precast_name] <- colData(spe)$layer

## Convert from character to a factor
colData(spe)[nnSVG_precast_name] <- as.factor(colData(spe)[,nnSVG_precast_name])

#avoid limma make contrasts syntax error if using 1-9 labels
#colData(spe)[nnSVG_precast_name] <- paste0("clust", colData(spe)[,nnSVG_precast_name])

load(
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_pseudobulk", paste0(nnSVG_precast_name,".Rdata"))
)

registration_mod <-
    registration_model(spe_pseudo, covars = NULL)

block_cor <-
    registration_block_cor(spe_pseudo, registration_model = registration_mod)

results_enrichment <-
    registration_stats_enrichment(
        spe_pseudo,
        block_cor = block_cor,
        covars = NULL,
        gene_ensembl = "gene_id",
        gene_name = "gene_name"
    )

modeling_results <- list(
    "enrichment" = results_enrichment
)

###save modeling results list
save(
    modeling_results,
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_DE", paste0(nnSVG_precast_name,".Rdata"))
)

colData(spe)$spatialLIBD <- colData(spe)$registration_variable

sig_genes <- sig_genes_extract(
    n = 50,
    modeling_results = modeling_results,
    model_type = "enrichment",
    sce_layer = spe_pseudo
)

write.csv(sig_genes, file = here::here("processed-data", "11_differential_expression","pseudobulk", "nnSVG_precast_DE",
                                       paste0(nnSVG_precast_name, "_sig_genes_50.csv")), row.names = FALSE)

indices <- c()

indices <- append(indices, which(sig_genes$gene == "PVALB")) #DLPFC nat neuro previous marker for L4
indices <- append(indices, which(sig_genes$gene == "FABP7")) #DLPFC nat neuro previous marker for L1
indices <- append(indices, which(sig_genes$gene == "CCK")) #DLPFC nat neuro previous marker for L6
indices <- append(indices, which(sig_genes$gene == "KRT17")) #DLPFC nat neuro new marker for L6
indices <- append(indices, which(sig_genes$gene == "AQP4")) #DLPFC nat neuro new marker for L1
indices <- append(indices, which(sig_genes$gene == "HPCAL1")) #DLPFC nat neuro new marker for L2
indices <- append(indices, which(sig_genes$gene == "TRABD2A")) #DLPFC nat neuro new marker & known L5 marker
indices <- append(indices, which(sig_genes$gene == "RELN")) #known L1 marker
indices <- append(indices, which(sig_genes$gene == "FREM3")) #known L3 marker

pdf(file = here::here("plots", "11_differential_expression","pseudobulk", "nnSVG_precast_DE",
                      paste0("boxplots_", nnSVG_precast_name, ".pdf")),
    width = 8.5, height = 8)

for (i in indices) {
    layer_boxplot(
        i,
        sig_genes = sig_genes,
        short_title = TRUE,
        sce_layer = spe_pseudo,
        col_bkg_box = "grey80",
        col_bkg_point = "grey40",
        col_low_box = "violet",
        col_low_point = "darkviolet",
        col_high_box = "skyblue",
        col_high_point = "dodgerblue4",
        cex = 2,
        group_var = "layer",
        assayname = "logcounts"
    )
}

dev.off()

#volcano plots
thresh_fdr <- 0.05
thresh_logfc <- log2(1.5)
fdrs_gene_ids <- rowData(spe_pseudo)$gene_id
fdrs_gene_names <- rowData(spe_pseudo)$gene_name

df_list <- list()

pdf(file = here::here("plots", "11_differential_expression","pseudobulk", "nnSVG_precast_DE",
                      paste0("volcano_", nnSVG_precast_name, ".pdf")),
    width = 8.5, height = 8)

for (i in unique(colData(spe_pseudo)[["layer"]])) {
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
                    title = "nnSVG PRECAST dACC",
                    subtitle = paste0(i, " vs. all others")
    ))
}

dev.off()
