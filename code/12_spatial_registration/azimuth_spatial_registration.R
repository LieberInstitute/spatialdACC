setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("here")
library("spatialLIBD")
library("SingleCellExperiment")
library("stringr")

# get reference layer enrichment statistics
#k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
k = 9
nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", k)
load(file = here("processed-data", "11_differential_expression", "pseudobulk",
                 "nnSVG_precast_DE", paste0(nnSVG_precast_name,".Rdata")))

# load cell type info from azimuth
load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))
table(sce$Sample, sce$cellType_azimuth)

sce_modeling_results <- registration_wrapper(
    sce = sce,
    var_registration = "cellType_azimuth",
    var_sample_id = "Sample",
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
)

## extract t-statics and rename
registration_t_stats <- sce_modeling_results$enrichment[, grep("^t_stat", colnames(sce_modeling_results$enrichment))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))

dim(registration_t_stats)

## check out table
registration_t_stats[1:5, 1:5]

cor_layer <- layer_stat_cor(
    stats = registration_t_stats,
    modeling_results = modeling_results,
    model_type = "enrichment",
    top_n = 100
)

cor_layer
save(cor_layer, file = here("processed-data", "12_spatial_registration",paste0("azimuth",nnSVG_precast_name,".rds")))

cor_layer <- cor_layer[c(1,19,9,18,7,10,11,8,16,12,17,13,15,3,5,14,2,4,6),]
t_cor_layer <- t(cor_layer)

layer_colors <- c(
    "L2" = "#377EB8",
    "L3" = "#4DAF4A",
    "L5" = "#FFD700",
    "L6b" = "#c46200",
    "L6a" = "#FFC18A",
    "WM" = "#1A1A1A",
    "L1" = "#F0027F"
)

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "celltype_colors.Rdata"))

pdf(file = here::here("plots", "12_spatial_registration", "azimuth",
                      paste0("azimuth_",nnSVG_precast_name,"_heatmap.pdf")), width = 6, height = 4)
layer_stat_cor_plot(t_cor_layer, color_max = max(cor_layer),
                    reference_colors = celltype_colors,
                    query_colors = layer_colors,
                    cluster_rows = FALSE,
                    cluster_columns = FALSE)
dev.off()

anno <- annotate_registered_clusters(
    cor_stats_layer = cor_layer,
    confidence_threshold = 0.25,
    cutoff_merge_ratio = 0.25
)

anno
