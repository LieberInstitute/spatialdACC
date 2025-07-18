setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(scater)
library(RcppML)
library(ggspavis)
library(here)
#library(scRNAseq)
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

p_list <- vis_grid_clus(spe, "classification",
              spatial = F, ncol = 4,
              return_plots = T)
png(here("plots", "13_NMF", "spot_plots_NMF38_NMF61.png"), height=20, width=20, unit="in",res=300)
cowplot::plot_grid(plotlist = p_list, ncol = 4)
dev.off()

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

save(modeling_results, file=here("processed-data", "13_NMF", "DE_NMF38_61.Rdata"))

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

# create violin plots of DRD5 expression in NMF_38 and NMF_61
# each spot is a data point

spe$counts_DRD5 <- counts(spe)[which(rowData(spe)$gene_name=="DRD5"),]
spe$logcounts_DRD5 <- logcounts(spe)[which(rowData(spe)$gene_name=="DRD5"),]

pdf(file = here::here("plots", "13_NMF", "DRD5_violin_plots_NMF_factors.pdf"),
    width = 10, height = 10)

plotColData(spe, x = "classification", y = "counts_DRD5") +
    ggtitle("DRD5 Counts by Factor") +
    facet_wrap(~ spe$classification, scales = "free_x", nrow = 1)

plotColData(spe, x = "classification", y = "logcounts_DRD5") +
    ggtitle("DRD5 Logcounts by Factor") +
    facet_wrap(~ spe$classification, scales = "free_x", nrow = 1)

dev.off()

# create volcano plot of pairwise comparison of NMF_38 and NMF_61
#volcano plots
thresh_fdr <- 0.05
thresh_logfc <- log2(1.5)
fdrs_gene_ids <- rowData(spe_pseudo)$gene_id
fdrs_gene_names <- rowData(spe_pseudo)$gene_name

fdrs <- modeling_results[["pairwise"]][,paste0("fdr_", "NMF_38-NMF_61")]
logfc <- modeling_results[["pairwise"]][,paste0("logFC_", "NMF_38-NMF_61")]

sig <- (fdrs < thresh_fdr) & (abs(logfc) > thresh_logfc)
print(table(sig))

df_list <- data.frame(
    gene_name = modeling_results[["pairwise"]]$gene,
    logFC = -logfc,
    FDR = fdrs,
    sig = sig
    )



# loading snRNAseq data to get colors

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))

#plot harmony UMAP colored by azimuth labels

# remove Sst Chodl
sce <- sce[,-which(sce$cellType_azimuth == "Sst Chodl")]
sce$cellType_azimuth <- as.factor(sce$cellType_azimuth)

mycolors <- pals::cols25()[1:19]
celltype_colors <- setNames(mycolors, unique(sce$cellType_azimuth))
celltype_colors  
#  L2_3_IT      Oligo        L6b      L5_IT      Pvalb        OPC      Lamp5 
#  "#1F78C8"  "#ff0000"  "#33a02c"  "#6A33C2"  "#ff7f00"  "#565656"  "#FFD700" 
#        Vip      L6_CT   MicroPVM        Sst      L6_IT       VLMC      Astro 
#  "#a6cee3"  "#FB6496"  "#b2df8a"  "#CAB2D6"  "#FDBF6F"  "#999999"  "#EEE685" 
#      L5_ET       Sncg L6_IT_Car3    L5_6_NP       Endo 
#  "#C8308C"  "#FF83FA"  "#C814FA"  "#0000FF"  "#36648B" 

# NMF61 = L5_ET
# NMF38 = L5_IT

# create a color vector for the volcano plot
keyvals <- ifelse(
    df_list$logFC < -1.5, celltype_colors["L5_IT"], # NMF38
    ifelse(df_list$logFC > 1.5, celltype_colors["L5_ET"], # NMF61
           'grey30')) # not significant
keyvals[is.na(keyvals)] <- 'grey30'
names(keyvals)[keyvals == celltype_colors["L5_IT"]] <- 'NMF38'
names(keyvals)[keyvals == celltype_colors["L5_ET"]] <- 'NMF61'
names(keyvals)[keyvals == 'grey30'] <- 'not sig.'


png(file = here::here("plots", "13_NMF", "volcano_plots_NMF38_NMF61.png"),
    width = 4, height = 4.5, res = 300, unit = "in")

EnhancedVolcano(df_list,
                lab = df_list$gene_name,
                x = 'logFC',
                y = 'FDR',
                xlim = c(-3, 5),
                ylim = c(0, -log10(10e-22)),
                legendPosition = "bottom",
                selectLab = c("VAT1L", "POU3F1", "SULF2", "FEZF2", "GABRQ"),
                FCcutoff = 1.5,
                pCutoff = 0.05,
                labSize = 4.5,
                ylab = "-log10 FDR",
                legendLabels = c('Not sig.','LogFC','FDR',
                                'FDR & LogFC'),
                title = "NMF38 < ----- > NMF61",
                subtitle = "",
                caption = "",
                drawConnectors=TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                colCustom=keyvals
                ) +
ggpubr::theme_pubr() +
theme(legend.position = "none") + guides(fill="none")
    

dev.off()
