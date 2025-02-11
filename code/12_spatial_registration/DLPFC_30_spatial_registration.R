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
spe <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium")
# create spatial labels for DLPFC_30
spe$BayesSpace_harmony_09[spe$BayesSpace_harmony_09 == 3] <- "L2"
spe$BayesSpace_harmony_09[spe$BayesSpace_harmony_09 == 8] <- "L4"
spe$BayesSpace_harmony_09[spe$BayesSpace_harmony_09 == 7] <- "L6"
spe$BayesSpace_harmony_09[spe$BayesSpace_harmony_09 == 5] <- "L3"
spe$BayesSpace_harmony_09[spe$BayesSpace_harmony_09 == 6] <- "WM"
spe$BayesSpace_harmony_09[spe$BayesSpace_harmony_09 == 4] <- "L5"
spe$BayesSpace_harmony_09[spe$BayesSpace_harmony_09 == 2] <- "L1"
spe$BayesSpace_harmony_09[spe$BayesSpace_harmony_09 == 1] <- "meninges"
spe$BayesSpace_harmony_09[spe$BayesSpace_harmony_09 == 9] <- "WM"
# remove meninges
spe <- spe[ , which(spe$BayesSpace_harmony_09 != "meninges")]

table(spe$sample_id, spe$BayesSpace_harmony_09)

spe_pseudo <-
    registration_pseudobulk(spe,
                            var_registration = "BayesSpace_harmony_09",
                            var_sample_id = "sample_id"
    )

spe_modeling_results <- registration_wrapper(
    spe,
    var_registration = "BayesSpace_harmony_09",
    var_sample_id = "sample_id",
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
)

save(spe_modeling_results, file = here("processed-data", "12_spatial_registration",paste0("DLPFC_30_DE",".Rdata")))

sig_genes <- sig_genes_extract(
    n = 30,
    modeling_results = spe_modeling_results,
    model_type = "enrichment",
    sce_layer = spe_pseudo
)

write.csv(sig_genes, file = here::here("processed-data", "11_differential_expression","pseudobulk", "nnSVG_precast_DE",
                                       paste0("DLPFC_30", "_sig_genes_30.csv")), row.names = FALSE)


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
save(cor_layer, file = here("processed-data", "12_spatial_registration",paste0("DLPFC_30_",nnSVG_precast_name,".rds")))

cor_layer <- cor_layer[c(6,2,3,1,4,5,7),]

layer_colors <- c(
    "L2" = "#377EB8",
    "L3" = "#4DAF4A",
    "L5" = "#FFD700",
    "L6b" = "#c46200",
    "L6a" = "#FFC18A",
    "WM" = "#1A1A1A",
    "L1" = "#F0027F"
)

ref_colors <- c(
    "L2" = "#377EB8",
    "L3" = "#4DAF4A",
    "L5" = "#FFD700",
    "L4" = "#984EA3",
    "L6" = "#FF7F00",
    "WM" = "#1A1A1A",
    "L1" = "#F0027F"
)

pdf(file = here::here("plots", "12_spatial_registration", "DLPFC_manual",
                      paste0("DLPFC_30_",nnSVG_precast_name,"_heatmap.pdf")), width = 5, height = 4)
layer_stat_cor_plot(t(cor_layer), color_max = max(cor_layer),
                    query_colors = layer_colors,
                    reference_colors = ref_colors,
                    cluster_rows = FALSE,
                    cluster_columns = FALSE)
dev.off()

anno <- annotate_registered_clusters(
    cor_stats_layer = cor_layer,
    confidence_threshold = 0.25,
    cutoff_merge_ratio = 0.25
)

anno
