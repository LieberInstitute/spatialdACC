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
library(EnhancedVolcano)

# load spe with NMF results
load(file = here("processed-data", "13_NMF", "spe_NMF.Rdata"))

# use "nmf_35" and "nmf_15" as the column names
# if the entry is larger than 0, classify as "NMF_35" or "NMF_15" respectively
# we want to create one new column - either "NMF_35" or "NMF_15" or "neither"
# we don't want the spots that overlap to be classified as both NMF_35 and NMF_15, these should be "neither" too
colData(spe)$classification <- with(colData(spe),
                                    ifelse(nmf35 > 0 & nmf15 > 0, "neither",         # both greater than 0
                                           ifelse(nmf35 > 0, "NMF_35",                 # only nmf_35 greater than 0
                                                  ifelse(nmf15 > 0, "NMF_15", "neither") # only nmf_15 greater than 0
                                           )
                                    )
)

table(colData(spe)$classification)
# neither  NMF_15  NMF_35
# 71148    1622     597

# remove spots that are classified as "neither" from the spe
spe <- spe[, colData(spe)$classification != "neither"]

table(colData(spe)$classification, colData(spe)$layer)

vis_grid_clus(spe, "classification",
              pdf_file = here("plots", "13_NMF", "spot_plots_NMF35_NMF15.pdf"),
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
                                        "NMF35_NMF15_sig_genes_50.csv"), row.names = FALSE)
write.csv(sig_genes_reverse, file = here::here("processed-data", "13_NMF",
                                               "NMF35_NMF15_sig_genes_reverse_50.csv"), row.names = FALSE)


# create volcano plot of pairwise comparison of NMF_35 and NMF_15
#volcano plots
thresh_fdr <- 0.05
thresh_logfc <- log2(1.5)
fdrs_gene_ids <- rowData(spe_pseudo)$gene_id
fdrs_gene_names <- rowData(spe_pseudo)$gene_name

fdrs <- modeling_results[["pairwise"]][,paste0("fdr_", "NMF_15-NMF_35")]
logfc <- modeling_results[["pairwise"]][,paste0("logFC_", "NMF_15-NMF_35")]

sig <- (fdrs < thresh_fdr) & (abs(logfc) > thresh_logfc)
print(table(sig))

df_list <- data.frame(
    gene_name = modeling_results[["pairwise"]]$gene,
    logFC = logfc,
    FDR = fdrs,
    sig = sig
    )

pdf(file = here::here("plots", "13_NMF", "volcano_plots_NMF15_NMF35.pdf"),
    width = 8.5, height = 8)

print(EnhancedVolcano(df_list,
                      lab = df_list$gene_name,
                      x = 'logFC',
                      y = 'FDR',
                      FCcutoff = 1.5,
                      pCutoff = 0.05,
                      ylab = "-log10 FDR",
                      legendLabels = c('Not sig.','Log (base 2) FC','FDR',
                                       'FDR & Log (base 2) FC'),
                      title = "nnSVG PRECAST dACC",
                      subtitle = "NMF15 minus NMF35",
                      )
      )

dev.off()






