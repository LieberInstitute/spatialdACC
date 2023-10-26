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

read_barcoded_csv <- function(x) {
    df <- read.csv(x)
    colnames(df) <- tolower(colnames(df))

    if (colnames(df)[2] == "cluster") {
        colnames(df)[2] <-
            gsub("gene_expression_", "", basename(dirname(x)))
    }
    return(df)
}

nnSVG_PRECAST_import <- function(spe, cluster_dir = file.path(tempdir(), "exported_clusters")) {
    clustering_files <-
        list.files(
            here::here("processed-data", "08_clustering", "PRECAST", nnSVG_precast_name),
            pattern = "clusters.csv",
            all.files = TRUE,
            full.names = TRUE,
            recursive = TRUE
        )
    clusters_list <- lapply(clustering_files, read_barcoded_csv)
    clusters <- Reduce(function(...) merge(..., by = "key", all = TRUE), clusters_list)
    cluster_cols <- which(colnames(clusters) != "key")
    colnames(clusters)[cluster_cols] <- paste0("", colnames(clusters)[cluster_cols])

    colData(spe)[,nnSVG_precast_name] <- clusters[,2]
    return(spe)
}

load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))

k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#import precast nnSVG clusters
#i re wrote the key when using PRECAST, just add the cluster columns manually
nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", k)
spe <- nnSVG_PRECAST_import(
    spe,
    cluster_dir = here::here("processed-data", "08_clustering", "PRECAST", nnSVG_precast_name)
)

## Convert from character to a factor
colData(spe)[nnSVG_precast_name] <- as.factor(colData(spe)[,nnSVG_precast_name])

spe_pseudo <-
    registration_pseudobulk(spe,
                            var_registration = nnSVG_precast_name,
                            var_sample_id = "sample_id",
                            min_ncells = 10
    )

registration_mod <-
    registration_model(spe_pseudo,
                       var_registration = nnSVG_precast_name)

block_cor <-
    registration_block_cor(spe_pseudo, registration_model = registration_mod)

results_enrichment <-
    registration_stats_enrichment(
        spe_pseudo,
        block_cor = block_cor,
        gene_ensembl = 'gene_id',
        gene_name = 'gene_name'
    )

results_anova <-
    registration_stats_anova(
        spe_pseudo,
        block_cor = block_cor,
        var_registration = nnSVG_precast_name,
        gene_ensembl = 'gene_id',
        gene_name = 'gene_name'
    )

modeling_results <- list(
    "anova" = results_anova,
    "enrichment" = results_enrichment
)

###save modeling results list

save(
    modeling_results,
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_DE", paste0(nnSVG_precast_name,".Rdata"))
)

thresh_fdr <- 0.05
thresh_logfc <- log2(1.5)
fdrs_gene_ids <- rowData(spe_pseudo)$gene_id
fdrs_gene_names <- rowData(spe_pseudo)$gene_name

df_list <- list()

pdf(file = here::here("plots", "11_differential_expression","pseudobulk", "nnSVG_precast_DE",
                      paste0("volcano_", nnSVG_precast_name, ".pdf")),
    width = 8.5, height = 8)


for (i in c(1:k)) {

    print(i)
    fdrs <- results_enrichment[,paste0("fdr_", i)]
    logfc <- results_enrichment[,paste0("logFC_", i)]

    # Identify significant genes (low FDR and high logFC)
    sig <- (fdrs < thresh_fdr) & (abs(logfc) > thresh_logfc)

    # Number of significant genes
    print(paste("Cluster", i))
    print(table(sig))

    df_list[[i]] <- data.frame(
        gene_name = results_enrichment$gene,
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
                    subtitle = paste0("Cluster ", i, " vs. all others")
    ))
}


