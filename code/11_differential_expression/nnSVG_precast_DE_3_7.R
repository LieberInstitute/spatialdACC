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

#avoid limma make contrasts syntax error
colData(spe)[nnSVG_precast_name] <- paste0("clust", colData(spe)[,nnSVG_precast_name])

spe_pseudo <-
    registration_pseudobulk(spe,
                            var_registration = nnSVG_precast_name,
                            var_sample_id = "sample_id",
                            min_ncells = 10
    )

modeling_results <- registration_wrapper(
    spe,
    var_registration = nnSVG_precast_name,
    var_sample_id = "sample_id",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

###load modeling results list
load(
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_DE", paste0(nnSVG_precast_name,".Rdata"))
)

colData(spe)$spatialLIBD <- colData(spe)$registration_variable

sig_genes <- sig_genes_extract(
    n = 20,
    modeling_results = modeling_results,
    model_type = "pairwise",
    sce_layer = spe_pseudo
)

## subset to cluster 3 vs cluster 7 comparison
sig_genes <- sig_genes[sig_genes$test == "clust3-clust7",]

indices <- c(1:20)

pdf(file = here::here("plots", "11_differential_expression","pseudobulk", "nnSVG_precast_DE",
                      paste0("clust3_clust7_boxplots_", nnSVG_precast_name, ".pdf")),
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
        group_var = nnSVG_precast_name,
        assayname = "logcounts"
    )
}

dev.off()

# save csv file of sig_genes
write.csv(sig_genes, file = here::here("processed-data", "11_differential_expression","pseudobulk", "nnSVG_precast_DE",
                                       paste0("clust3_clust7_", nnSVG_precast_name, "_sig_genes.csv")), row.names = FALSE)

