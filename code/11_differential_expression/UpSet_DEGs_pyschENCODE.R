setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("scater")
    library("spatialLIBD")
    library("dplyr")
    library('UpSetR')
})

nnSVG_precast_name <- paste0("nnSVG_PRECAST_captureArea_", 9)
sig_genes_dACC <- read.csv(file = here::here("processed-data", "11_differential_expression","pseudobulk", "nnSVG_precast_DE",
                                             paste0(nnSVG_precast_name, "_sig_genes_30.csv")))

sig_genes_DLPFC <- read.csv(file = here::here("processed-data", "11_differential_expression","pseudobulk", "nnSVG_precast_DE",
                                              paste0("DLPFC_30", "_sig_genes_30.csv")))

# remove WM genes from both lists
sig_genes_dACC <- sig_genes_dACC[sig_genes_dACC$test != "WM",]
sig_genes_DLPFC <- sig_genes_DLPFC[sig_genes_DLPFC$test != "WM",]

layers_dACC <- unique(sig_genes_dACC$test)
layers_DLPFC <- unique(sig_genes_DLPFC$test)
layers_dACC; layers_DLPFC

gene_sets <- list()

# separate comparisons for L4 with L3 and L5
gene_sets[["L4_DLPFC"]] <- unique(sig_genes_DLPFC$gene[sig_genes_DLPFC$test == "L4"])
gene_sets[["L3_dACC"]] <- unique(sig_genes_dACC$gene[sig_genes_dACC$test == "L3"])
gene_sets[["L5_dACC"]] <- unique(sig_genes_dACC$gene[sig_genes_dACC$test == "L5"])

# separate comparisons for L6 with L6a and L6b
gene_sets[["L6_DLPFC"]] <- unique(sig_genes_DLPFC$gene[sig_genes_DLPFC$test == "L6"])
gene_sets[["L6a_dACC"]] <- unique(sig_genes_dACC$gene[sig_genes_dACC$test == "L6a"])
gene_sets[["L6b_dACC"]] <- unique(sig_genes_dACC$gene[sig_genes_dACC$test == "L6b"])

for (layer in intersect(layers_dACC, layers_DLPFC)) {
    if (!layer %in% c("L4", "L6")) {
        gene_sets[[paste0(layer, "_dACC")]] <- unique(sig_genes_dACC$gene[sig_genes_dACC$test == layer])
        gene_sets[[paste0(layer, "_DLPFC")]] <- unique(sig_genes_DLPFC$gene[sig_genes_DLPFC$test == layer])
    }
}

upset_data <- fromList(gene_sets)

pdf(file = here::here("plots", "11_differential_expression", "upset_dACC_DLPFC_30.pdf"), width = 10, height = 10)
upset(upset_data,
      sets = names(gene_sets),
      main.bar.color = "blue",
      matrix.color = "darkred",
      order.by = "freq",
      keep.order = TRUE)
dev.off()
