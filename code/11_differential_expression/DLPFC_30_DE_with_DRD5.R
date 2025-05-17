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

save(spe_pseudo, file = here("processed-data", "12_spatial_registration",paste0("DLPFC_30_pseudobulk",".Rdata")))
assay(spe_pseudo, "logcounts")<- NULL

load(file = here("processed-data", "12_spatial_registration",paste0("DLPFC_30_pseudobulk",".Rdata")))
spatialCoords(DRD5_DLPFC_30_pseudo) <- NULL

spe_pseudo <- rbind(spe_pseudo, DRD5_DLPFC_30_pseudo[6497,])

logcounts(spe_pseudo) <-
    edgeR::cpm(
        edgeR::calcNormFactors(spe_pseudo),
        log = TRUE,
        prior.count = 1
    )

if (is(spe_pseudo, "SpatialExperiment")) {
    spatialCoords(spe_pseudo) <- NULL
    imgData(spe_pseudo) <- NULL
}

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

results_enrichment <-
    registration_stats_enrichment(
        spe_pseudo,
        block_cor = block_cor,
        covars = NULL,
        gene_ensembl = "gene_id",
        gene_name = "gene_name"
    )

results_anova <-
    registration_stats_anova(
        spe_pseudo,
        block_cor = block_cor,
        covars = NULL,
        gene_ensembl = "gene_id",
        gene_name = "gene_name"
    )

modeling_results <- list(
    "pairwise" = results_pairwise,
    "enrichment" = results_enrichment,
    "anova" = results_anova
)

save(modeling_results, file = here("processed-data", "11_differential_expression",paste0("DLPFC_30_DE_with_DRD5",".Rdata")))
