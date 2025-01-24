setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("here")
library("spatialLIBD")
library("SingleCellExperiment")
library("stringr")

k=9
nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", k)

modeling_results_manual <- readRDS(file = here("processed-data", "20_WM_comparisons", "modeling_results_WM_vs_CC.RDS"))

## extract t-statics and rename
registration_t_stats <- modeling_results_manual$enrichment[, grep("^t_stat", colnames(modeling_results_manual$enrichment))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))

load(file = here("processed-data", "12_spatial_registration",paste0("DLPFC_30_DE",".Rdata")))

cor_layer <- layer_stat_cor(
    stats = registration_t_stats,
    modeling_results = spe_modeling_results,
    model_type = "enrichment",
    top_n = 100
)

cor_layer
save(cor_layer, file = here("processed-data", "12_spatial_registration",paste0("dACC_manual",nnSVG_precast_name,".rds")))

cor_layer <- cor_layer[c(6,7,8,4,5,3,2,1),]

pdf(file = here::here("plots", "12_spatial_registration", "dACC_manual",
                      paste0("dACC_manual_","DLPFC_30","_heatmap.pdf")), width = 5, height = 5)
layer_stat_cor_plot(cor_layer, max = max(cor_layer))
title("Manual vs. DLPFC")
dev.off()

anno <- annotate_registered_clusters(
    cor_stats_layer = cor_layer,
    confidence_threshold = 0.25,
    cutoff_merge_ratio = 0.25
)

anno
