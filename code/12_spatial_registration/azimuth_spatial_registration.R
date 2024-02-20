setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("here")
library("spatialLIBD")
library("SingleCellExperiment")
library("stringr")

# get reference layer enrichment statistics
k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", k)
load(file = here("processed-data", "11_differential_expression", "pseudobulk",
                 "nnSVG_precast_DE", paste0(nnSVG_precast_name,".Rdata")))

modeling_results$enrichment[1:5, 1:5]

# load cell type info from azimuth
load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))
table(sce$Sample, sce$cellType_azimuth)

# replace "-" with nothing
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "Micro-PVM", "MicroPVM")

# replace spaces with underscores
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L2/3 IT", "L2_3_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5 ET", "L5_ET")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5 IT", "L5_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5/6 NP", "L5_6_NP")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 CT", "L6_CT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 IT", "L6_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 IT Car3", "L6_IT_Car3")

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

## cell types x gene
dim(registration_t_stats)
#> [1] 21527    19

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

pdf(file = here::here("plots", "12_spatial_registration", "azimuth",
                      paste0("azimuth_",nnSVG_precast_name,"_heatmap.pdf")), width = 14, height = 14)
layer_stat_cor_plot(cor_layer, max = max(cor_layer))
dev.off()

anno <- annotate_registered_clusters(
    cor_stats_layer = cor_layer,
    confidence_threshold = 0.25,
    cutoff_merge_ratio = 0.25
)

anno
