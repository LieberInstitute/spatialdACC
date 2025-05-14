setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("here")
library("spatialLIBD")
library("SingleCellExperiment")
library("stringr")

k=9
nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", k)
load(file = here("processed-data", "11_differential_expression", "pseudobulk",
                 "nnSVG_precast_DE", paste0(nnSVG_precast_name,".Rdata")))

modeling_results_manual <- readRDS(file = here("processed-data", "20_WM_comparisons", "modeling_results_WM_vs_CC.RDS"))

## extract t-statics and rename
registration_t_stats <- modeling_results_manual$enrichment[, grep("^t_stat", colnames(modeling_results_manual$enrichment))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))

cor_layer <- layer_stat_cor(
    stats = registration_t_stats,
    modeling_results = modeling_results,
    model_type = "enrichment",
    top_n = 100
)

cor_layer
save(cor_layer, file = here("processed-data", "12_spatial_registration",paste0("dACC_manual",nnSVG_precast_name,".rds")))

cor_layer <- cor_layer[c(5,6,3,4,2,1),]

layer_colors <- c(
    "L2" = "#377EB8",
    "L3" = "#4DAF4A",
    "L5" = "#FFD700",
    "L6b" = "#c46200",
    "L6a" = "#FFC18A",
    "WM" = "#1A1A1A",
    "L1" = "#F0027F"
)

pdf(file = here::here("plots", "12_spatial_registration", "dACC_manual",
                      paste0("dACC_manual_",nnSVG_precast_name,"_heatmap.pdf")), width = 6, height = 4)
layer_stat_cor_plot(t(cor_layer), color_max = max(cor_layer),
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
