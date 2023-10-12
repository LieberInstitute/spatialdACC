setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("scater")
    library("spatialLIBD")
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
dim(spe_pseudo)

pca <- prcomp(t(assays(spe_pseudo)$logcounts))
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
pca_pseudo <- pca$x[, seq_len(20)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

## save pseudobulked spe file
save(
    spe_pseudo,
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_pseudobulk", paste0(nnSVG_precast_name,".Rdata"))
)

## Plot PCs
colData(spe_pseudo)["layer"] <- colData(spe_pseudo)[,nnSVG_precast_name]

pdf(file = here("plots", "11_differential_expression", "pseudobulk", "nnSVG_precast_pseudobulk", paste0("pseudobulk_PC_",nnSVG_precast_name,".pdf")), width = 14, height = 14)

plotPCA(
        spe_pseudo,
        colour_by = "layer",
        ncomponents = 2,
        point_size = 2,
        percentVar = metadata(spe_pseudo)$PCA_var_explained
    )

plotPCA(
    spe_pseudo,
    colour_by = "sample_id",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(spe_pseudo)$PCA_var_explained
)

plotPCA(
        spe_pseudo,
        colour_by = "sum",
        ncomponents = 2,
        point_size = 2,
        percentVar = metadata(spe_pseudo)$PCA_var_explained
    )

plotPCA(
    spe_pseudo,
    colour_by = "detected",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(spe_pseudo)$PCA_var_explained
)
vars <- getVarianceExplained(spe_pseudo,
                             variables = c("layer","sample_id")
)


plotExplanatoryVariables(vars)

dev.off()
