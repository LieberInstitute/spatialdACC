setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(scater)
library(RcppML)
library(ggspavis)
library(here)
library(scRNAseq)
library(Matrix)
library(scran)
library(scuttle)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(igraph)
library(bluster)
library(patchwork)
library(cowplot)
library(projectR)
library(spatialLIBD)
library(gridExtra)
library(escheR)

# load spe with NMF results
load(file = here("processed-data", "13_NMF", "spe_NMF.Rdata"))

# use "nmf_38" and "nmf_61" as the column names
# if the entry is larger than 0, classify as "NMF_38" or "NMF_61" respectively
# we want to create one new column - either "NMF_38" or "NMF_61" or "neither"
# we don't want the spots that overlap to be classified as both NMF_38 and NMF_61, these should be "neither" too
colData(spe)$classification <- with(colData(spe),
                                    ifelse(nmf38 > 0 & nmf61 > 0, "neither",         # both greater than 0
                                           ifelse(nmf38 > 0, "NMF_38",                 # only nmf_38 greater than 0
                                                  ifelse(nmf61 > 0, "NMF_61", "neither") # only nmf_61 greater than 0
                                           )
                                    )
)

table(colData(spe)$classification)
# neither  NMF_38  NMF_61
# 68154    4191    1022

# remove spots that are classified as "neither" from the spe
spe <- spe[, colData(spe)$classification != "neither"]

table(colData(spe)$classification, colData(spe)$layer)

vis_grid_clus(spe, "classification",
              pdf_file = here("plots", "13_NMF", "spot_plots_NMF38_NMF61.pdf"),
              spatial = F, ncol = 3)

spe_pseudo <-
    registration_pseudobulk(spe,
                            var_registration = "classification",
                            var_sample_id = "sample_id",
                            min_ncells = 10
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

modeling_results <- list(
    "pairwise" = results_pairwise
)

sig_genes <- sig_genes_extract(
    n = 50,
    modeling_results = modeling_results,
    model_type = "pairwise",
    sce_layer = spe_pseudo
)

sig_genes_reverse <- sig_genes_extract(
    n = 50,
    modeling_results = modeling_results,
    model_type = "pairwise",
    reverse = TRUE,
    sce_layer = spe_pseudo
)

write.csv(sig_genes, file = here::here("processed-data", "13_NMF",
                                        "NMF38_NMF61_sig_genes_50.csv"), row.names = FALSE)
write.csv(sig_genes_reverse, file = here::here("processed-data", "13_NMF",
                                               "NMF38_NMF61_sig_genes_reverse_50.csv"), row.names = FALSE)
