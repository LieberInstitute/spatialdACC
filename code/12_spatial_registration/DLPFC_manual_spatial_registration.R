setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("here")
library("spatialLIBD")
library("SingleCellExperiment")
library("stringr")

# get reference layer enrichment statistics
#k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

k=9
nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", k)
load(file = here("processed-data", "11_differential_expression", "pseudobulk",
                 "nnSVG_precast_DE", paste0(nnSVG_precast_name,".Rdata")))

# load DLPFC manual annotations
spe <- spatialLIBD::fetch_data(type = "spe")

spe$layer_guess_reordered <- unfactor(spe$layer_guess_reordered)
spe$layer_guess_reordered[spe$layer_guess_reordered == "Layer1"] <- "L1"
spe$layer_guess_reordered[spe$layer_guess_reordered == "Layer2"] <- "L2"
spe$layer_guess_reordered[spe$layer_guess_reordered == "Layer3"] <- "L3"
spe$layer_guess_reordered[spe$layer_guess_reordered == "Layer4"] <- "L4"
spe$layer_guess_reordered[spe$layer_guess_reordered == "Layer5"] <- "L5"
spe$layer_guess_reordered[spe$layer_guess_reordered == "Layer6"] <- "L6"
# remove NA
spe <- spe[,!is.na(spe$layer_guess_reordered)]

table(spe$sample_id, spe$layer_guess_reordered)

spe_pseudo <-
    registration_pseudobulk(spe,
                            var_registration = "layer_guess_reordered",
                            var_sample_id = "sample_id"
    )

spe_modeling_results <- registration_wrapper(
    spe,
    var_registration = "layer_guess_reordered",
    var_sample_id = "sample_id",
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
)

sig_genes <- sig_genes_extract(
    n = 30,
    modeling_results = spe_modeling_results,
    model_type = "enrichment",
    sce_layer = spe_pseudo
)

write.csv(sig_genes, file = here::here("processed-data", "11_differential_expression","pseudobulk", "nnSVG_precast_DE",
                                       paste0("DLPFC_12", "_sig_genes_30.csv")), row.names = FALSE)


## extract t-statics and rename
registration_t_stats <- spe_modeling_results$enrichment[, grep("^t_stat", colnames(spe_modeling_results$enrichment))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))

## cell types x gene
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
save(cor_layer, file = here("processed-data", "12_spatial_registration",paste0("DLPFC_12_",nnSVG_precast_name,".rds")))

pdf(file = here::here("plots", "12_spatial_registration", "DLPFC_manual",
                      paste0("DLPFC_12_",nnSVG_precast_name,"_heatmap.pdf")), width = 14, height = 14)
layer_stat_cor_plot(cor_layer, max = max(cor_layer))
title("DLPFC n=12 vs. dACC")
dev.off()

anno <- annotate_registered_clusters(
    cor_stats_layer = cor_layer,
    confidence_threshold = 0.25,
    cutoff_merge_ratio = 0.25
)

anno
