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
                                             paste0(nnSVG_precast_name, "_sig_genes_50.csv")))

sig_genes_DLPFC <- read.csv(file = here::here("processed-data", "11_differential_expression","pseudobulk", "nnSVG_precast_DE",
                                              paste0("DLPFC_30", "_sig_genes_50.csv")))

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

pdf(file = here::here("plots", "11_differential_expression", "upset_50_DEGs_dACC_DLPFC_30.pdf"), width = 10, height = 10)
upset(upset_data,
      sets = names(gene_sets),
      main.bar.color = "blue",
      matrix.color = "darkred",
      order.by = "freq",
      keep.order = TRUE)
dev.off()

upset_data_table <- data.frame(gene = character(), layer = character(), stringsAsFactors = FALSE)

add_combination_to_table <- function(gene_sets, genes, layers) {
    upset_data_table <<- rbind(upset_data_table, data.frame(gene = genes, layer = layers, stringsAsFactors = FALSE))
}

all_combinations <- combn(names(gene_sets), 2, simplify = FALSE)  # pairwise combinations

for (comb in all_combinations) {
    intersection_genes <- Reduce(intersect, gene_sets[comb])
    if (length(intersection_genes) > 0) {
        add_combination_to_table(gene_sets, intersection_genes, paste(comb, collapse = " & "))
    }
}

for (k in 3:length(gene_sets)) {
    all_combinations_k <- combn(names(gene_sets), k, simplify = FALSE)

    for (comb in all_combinations_k) {
        intersection_genes <- Reduce(intersect, gene_sets[comb])
        if (length(intersection_genes) > 0) {
            add_combination_to_table(gene_sets, intersection_genes, paste(comb, collapse = " & "))
        }
    }
}

added_genes <- unique(upset_data_table$gene)  # track genes already added

# add  non-overlapping genes for each layer
for (layer in names(gene_sets)) {
    non_overlap_genes <- setdiff(gene_sets[[layer]], added_genes)
    if (length(non_overlap_genes) > 0) {
        add_combination_to_table(gene_sets, non_overlap_genes, layer)
        added_genes <- c(added_genes, non_overlap_genes)
    }
}

duplicates <- upset_data_table[duplicated(upset_data_table$gene), ]
for (gene in unique(duplicates$gene)) {
    print(gene)
    indices_gene <- which(upset_data_table$gene == gene)
    if (length(indices_gene) > 1) {
        upset_data_table <- upset_data_table[-c(indices_gene[-length(indices_gene)]),]
    }
}

duplicates <- upset_data_table[duplicated(upset_data_table$gene), ]
# none left


write.csv(upset_data_table, here("processed-data", "11_differential_expression", "UpSet_table_50_DEGs_DLPFC_30.csv"), row.names = FALSE)
