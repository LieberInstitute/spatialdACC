setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("scater")
    library("spatialLIBD")
    library("dplyr")
    library('EnhancedVolcano')
    library("patchwork")
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

results_pairwise <-
    registration_stats_pairwise(
        spe_pseudo,
        registration_model = registration_mod,
        block_cor = block_cor,
        gene_ensembl = "gene_id",
        gene_name = "gene_name"
    )

results_enrichment <-
    registration_stats_enrichment(
        spe_pseudo,
        block_cor = block_cor,
        covars = NULL,
        gene_ensembl = "gene_id",
        gene_name = "gene_name"
    )

results_anova <-
    registration_stats_anova(
        spe_pseudo,
        block_cor = block_cor,
        covars = NULL,
        gene_ensembl = "gene_id",
        gene_name = "gene_name"
    )

modeling_results <- list(
    "pairwise" = results_pairwise,
    "enrichment" = results_enrichment,
    "anova" = results_anova
)

###save modeling results list
save(
    modeling_results,
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_DE", paste0(nnSVG_precast_name,".Rdata"))
)

colData(spe)$spatialLIBD <- colData(spe)$registration_variable

sig_genes <- sig_genes_extract(
    n = 30,
    modeling_results = modeling_results,
    model_type = "enrichment",
    sce_layer = spe_pseudo
)

write.csv(sig_genes, file = here::here("processed-data", "11_differential_expression","pseudobulk", "nnSVG_precast_DE",
                                       paste0(nnSVG_precast_name, "_sig_genes_30.csv")), row.names = FALSE)


## For sig_genes_extract_all() to work
spe_pseudo$spatialLIBD <- spe_pseudo$layer

sig_genes_all <- sig_genes_extract_all(
    n = 30,
    modeling_results = modeling_results,
    sce_layer = spe_pseudo
)

saveRDS(sig_genes_all, file = here::here("processed-data", "11_differential_expression","pseudobulk", "nnSVG_precast_DE",
                                       paste0(nnSVG_precast_name, "_sig_genes_all.rds")))


# remove the "WM" rows in the "test" column
sig_genes <- sig_genes[sig_genes$test != "WM",]

# find genes that are repeated in the "gene" column
repeated_genes <- sig_genes[duplicated(sig_genes$gene),]
# GRIA4 in L2 and L3 & KRT17 in L6a and L6b
which(sig_genes$gene == "GRIA4")
which(sig_genes$gene == "KRT17")

indices <- c()

indices <- append(indices, which(sig_genes$gene == "GRIA4")) #new marker for L2 and L3
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
plot_list <- list()

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


# create volcano plot of pairwise comparison of L6a and L6b
#volcano plots
thresh_fdr <- 0.05
thresh_logfc <- log2(1.5)
fdrs_gene_ids <- rowData(spe_pseudo)$gene_id
fdrs_gene_names <- rowData(spe_pseudo)$gene_name

fdrs <- modeling_results[["pairwise"]][,paste0("fdr_", "L6a-L6b")]
logfc <- modeling_results[["pairwise"]][,paste0("logFC_", "L6a-L6b")]

sig <- (fdrs < thresh_fdr) & (abs(logfc) > thresh_logfc)
print(table(sig))

df_list <- data.frame(
    gene_name = modeling_results[["pairwise"]]$gene,
    logFC = logfc,
    FDR = fdrs,
    sig = sig
)

p_pair <- EnhancedVolcano(df_list,
                          lab = df_list$gene_name,
                          x = 'logFC',
                          y = 'FDR',
                          xlim = c(-3, 4),
                          ylim = c(0, -log10(10e-40)),
                          legendPosition = "bottom",
                          selectLab = c("CCK", "NPTXR", "NCAM2", "MBP", "MOG", "GFAP"),
                          FCcutoff = 1,
                          pCutoff = 0.05,
                          labSize = 7.0,
                          ylab = "-log10 FDR",
                          legendLabels = c('Not sig.','LogFC','FDR',
                                           'FDR & LogFC'),
                          title = "L6a vs. L6b",
                          subtitle = "",
                          caption = ""
)


# supp figure
png(file = here::here("plots", "11_differential_expression","pseudobulk", "nnSVG_precast_DE",
                      paste0("volcano_", nnSVG_precast_name, ".png")),
    width = 13, height = 18, unit="in", res=300)

wrap_plots(plot_list[["L1"]],plot_list[["L2"]],plot_list[["L3"]],plot_list[["L5"]],
           plot_list[["L6a"]],plot_list[["L6b"]],plot_list[["WM"]], p_pair,
           nrow=4, guides="collect") + plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom')

dev.off()

pdf(file = here::here("plots", "11_differential_expression","pseudobulk", "nnSVG_precast_DE",
                      paste0("volcano_", "L6a-L6b", ".pdf")),
    width = 7, height = 7)

print(p_pair)

dev.off()

sig_genes <- sig_genes_extract(
    n = 50,
    modeling_results = modeling_results,
    model_type = "pairwise",
    sce_layer = spe_pseudo
)

idx <- which(sig_genes$test == "L6a-L6b")

sig_genes <- sig_genes[idx,]

sig_genes_reverse <- sig_genes_extract(
    n = 50,
    modeling_results = modeling_results,
    model_type = "pairwise",
    reverse = TRUE,
    sce_layer = spe_pseudo
)

idx <- which(sig_genes_reverse$test == "L6b-L6a")

sig_genes_reverse <- sig_genes_reverse[idx,]

sig_genes_combined <- rbind(sig_genes,sig_genes_reverse)

write.csv(sig_genes_combined, file = here::here("processed-data", "11_differential_expression",
                                       "L6a_L6b_sig_genes_50.csv"), row.names = FALSE)
