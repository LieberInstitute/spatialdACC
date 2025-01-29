setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library(here)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(SingleCellExperiment)

# load DEGs
load(
    file = here("processed-data", "snRNA-seq", "05_azimuth", "azimuth_DE_results.Rdata")
)

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))

# replace spaces with underscores
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L2/3 IT", "L2_3_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5 ET", "L5_ET")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5 IT", "L5_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L5/6 NP", "L5_6_NP")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 CT", "L6_CT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 IT", "L6_IT")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "L6 IT Car3", "L6_IT_Car3")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "Sst Chodl", "Sst_Chodl")
colData(sce)$cellType_azimuth <- replace(colData(sce)$cellType_azimuth, colData(sce)$cellType_azimuth == "MicroPVM", "Micro_PVM")

# Subset modeling results
enrichment_results <- modeling_results[["enrichment"]]
dim(enrichment_results)
# 21527    78

# remove duplicates
enrichment_results <- enrichment_results[!duplicated(enrichment_results$gene),]

# Load the bulk data
file = here("processed-data", "pain", "GSE227159_Processed-data_Becker-Lutz-Yalcin_BLA-ACC-pathway.csv")
bulk <- read.csv(file, row.names = 2)

dim(bulk)
# [1] 37786    48

# better formatting
bulk$logFC <- bulk$Log2.FC..ctrl.vs.stim.
bulk$adj_p <- bulk$Adjusted.p.value..ctrl.vs.stim.
bulk$p <- bulk$P.value..ctrl.vs.stim.

bulk <- bulk[!(is.na(bulk$adj_p) | bulk$adj_p=="" | bulk$p==""), ]

bulk$logFC <- str_replace(bulk$logFC, ",", ".")
bulk$adj_p <- str_replace(bulk$adj_p, ",", ".")
bulk$p <- str_replace(bulk$p, ",", ".")

# List columns of interest from bulk data
vars <- c("Gene.name", "logFC", "adj_p", "p")

# Create smaller df with these vars
bulk <- bulk[, vars]

bulk$p <- formatC(bulk$p, format = "f", digits = 20)
bulk$p <- as.numeric(bulk$p)

orthology<-read.csv(file=here::here('raw-data','Retro-seq',
                                    'human_mouse_orthologs.csv'))

names <- orthology[orthology$Column3 %in% bulk$Gene.name,]

names <- names[match(bulk$Gene.name, names$Column3),]

setdiff(names$Column3, bulk$Gene.name)

sum(is.na(names$Column1))
# [1] 1436

bulk$Human.gene.name <- names$Column1

# remove non-matching human genes
bulk <- bulk[!is.na(bulk$Human.gene.name),]
dim(bulk)
# 12820 4

length(unique(bulk$Human.gene.name))
# [1] 12659

#bulk[duplicated(bulk$Human.gene.name),]
# these genes are not of interest since they are not significant DEGs,
# just remove one instance of each to avoid errors downstream

bulk <- bulk[!duplicated(bulk$Human.gene.name),]

# Overlap to get genes in both datasets
overlap <- intersect(bulk$Human.gene.name, enrichment_results$gene)
bulk <- bulk[bulk$Human.gene.name %in% overlap, ]
enrichment_results <- enrichment_results[enrichment_results$gene %in% overlap, ]

adj_pval_threshold <- 0.1

DE_bulk_up <- bulk[bulk$adj_p < adj_pval_threshold & bulk$logFC > 0, ]$Human.gene.name
DE_bulk_down <- bulk[bulk$adj_p < adj_pval_threshold & bulk$logFC < 0, ]$Human.gene.name

nonDE_bulk_up <- setdiff(bulk$Human.gene.name, DE_bulk_up)
nonDE_bulk_down <- setdiff(bulk$Human.gene.name, DE_bulk_down)

#pval_threshold <- 0.05

#DE_bulk_up <- bulk[bulk$p < pval_threshold & bulk$logFC > 0, ]$Human.gene.name
#DE_bulk_down <- bulk[bulk$p < pval_threshold & bulk$logFC < 0, ]$Human.gene.name

#nonDE_bulk_up <- setdiff(bulk$Human.gene.name, DE_bulk_up)
#nonDE_bulk_down <- setdiff(bulk$Human.gene.name, DE_bulk_down)


# Create a list to store the results
results_up <- list()
results_down <- list()

k <- unique(sce$cellType_azimuth)
k <- k[k != "Sst_Chodl"]

#top_n <- 100

for (i in k) {
    print(i)
    # Identify upregulated and downregulated DE genes in the spatial domain

    #t_stat_threshold <- sort(enrichment_results[[paste0("t_stat_", i)]], decreasing = T)[top_n]

    #DE_clust_genes_up <- enrichment_results[enrichment_results[[paste0("t_stat_", i)]] >= t_stat_threshold, ]$gene
    #DE_clust_genes_down <- DE_clust_genes_up

    #nonDE_clust_genes_up <- enrichment_results[enrichment_results[[paste0("t_stat_", i)]] < t_stat_threshold, ]$gene
    #nonDE_clust_genes_down <- nonDE_clust_genes_up

    DE_clust_genes_up <- enrichment_results[
        enrichment_results[[paste0("fdr_", i)]] < 0.05 & enrichment_results[[paste0("logFC_", i)]] > 1, ]$gene
    DE_clust_genes_down <- DE_clust_genes_up

    print(length(DE_clust_genes_up))

    nonDE_clust_genes_up <- enrichment_results[
        enrichment_results[[paste0("fdr_", i)]] >= 0.05 | enrichment_results[[paste0("logFC_", i)]] <= 1, ]$gene
    nonDE_clust_genes_down <- nonDE_clust_genes_up

    # Count overlaps for upregulated genes
    DE_clust_DE_bulk_up <- length(intersect(DE_clust_genes_up, DE_bulk_up))
    DE_clust_nonDE_bulk_up <- length(intersect(DE_clust_genes_up, nonDE_bulk_up))
    nonDE_clust_DE_bulk_up <- length(intersect(nonDE_clust_genes_up, DE_bulk_up))
    nonDE_clust_nonDE_bulk_up <- length(intersect(nonDE_clust_genes_up, nonDE_bulk_up))

    # Create the 2x2 table for upregulated genes
    contingency_table_up <- matrix(c(DE_clust_DE_bulk_up, DE_clust_nonDE_bulk_up,
                                     nonDE_clust_DE_bulk_up, nonDE_clust_nonDE_bulk_up),
                                   nrow = 2, byrow = TRUE,
                                   dimnames = list("Cluster" = c("DE", "nonDE"),
                                                   "Bulk" = c("DE", "nonDE")))

    # Perform the chi-square test and Fisher's exact test for upregulated genes
    chi_sq_test_up <- chisq.test(contingency_table_up)
    fisher_test_up <- fisher.test(contingency_table_up)

    # Store the results for upregulated genes
    results_up[[i]] <- list(contingency_table = contingency_table_up,
                            chi_sq_test = chi_sq_test_up,
                            fisher_test = fisher_test_up)

    # Count overlaps for downregulated genes
    DE_clust_DE_bulk_down <- length(intersect(DE_clust_genes_down, DE_bulk_down))
    DE_clust_nonDE_bulk_down <- length(intersect(DE_clust_genes_down, nonDE_bulk_down))
    nonDE_clust_DE_bulk_down <- length(intersect(nonDE_clust_genes_down, DE_bulk_down))
    nonDE_clust_nonDE_bulk_down <- length(intersect(nonDE_clust_genes_down, nonDE_bulk_down))

    # Create the 2x2 table for downregulated genes
    contingency_table_down <- matrix(c(DE_clust_DE_bulk_down, DE_clust_nonDE_bulk_down,
                                       nonDE_clust_DE_bulk_down, nonDE_clust_nonDE_bulk_down),
                                     nrow = 2, byrow = TRUE,
                                     dimnames = list("Cluster" = c("DE", "nonDE"),
                                                     "Bulk" = c("DE", "nonDE")))

    # Perform the chi-square test and Fisher's exact test for downregulated genes
    chi_sq_test_down <- chisq.test(contingency_table_down)
    fisher_test_down <- fisher.test(contingency_table_down)

    # Store the results for downregulated genes
    results_down[[i]] <- list(contingency_table = contingency_table_down,
                              chi_sq_test = chi_sq_test_down,
                              fisher_test = fisher_test_down)
}

# Print results for upregulated genes
for (i in k) {
    print(paste("Cluster", i, "- Upregulated Genes:"))
    print(results_up[[i]]$fisher_test)
}

# Print results for downregulated genes
for (i in k) {
    print(paste("Cluster", i, "- Downregulated Genes:"))
    print(results_down[[i]]$fisher_test)
}

# Initialize matrices for upregulated and downregulated p-values
pvalues_up <- matrix(NA, nrow = length(k), ncol = 1)
pvalues_down <- matrix(NA, nrow = length(k), ncol = 1)

# Fill matrices with p-values from Fisher's test
for (i in 1:length(k)) {
    index <- k[i]
    pvalues_up[i, ] <- results_up[[index]]$fisher_test$p.value
    pvalues_down[i, ] <- results_down[[index]]$fisher_test$p.value
}

# Combine matrices into one matrix for the heatmap
combined_pvalues <- cbind(
    Upreg = pvalues_up,
    Downreg = pvalues_down
)

# Assign row names for spatial domains and column names for regulation types
rownames(combined_pvalues) <- k
colnames(combined_pvalues) <- c("Upreg", "Downreg")

# Compute average values for reordering
row_means <- rowMeans(-log10(combined_pvalues), na.rm = TRUE)
col_means <- colMeans(-log10(combined_pvalues), na.rm = TRUE)

# Order rows and columns based on average values
ordered_rows <- order(row_means, decreasing = TRUE)
ordered_cols <- order(col_means, decreasing = TRUE)

# Reorder the heatmap matrix
combined_pvalues_ordered <- combined_pvalues[ordered_rows, ordered_cols]

col_fun <- colorRamp2(
    c(1.3, max(-log10(combined_pvalues_ordered))),
    c("white", "blue") # White for -log10(p) >= 1.3 (p >= 0.05), blue for more significant p-values
)

# Create heatmap for the combined p-values
heatmap_combined <- Heatmap(
    -log10(combined_pvalues_ordered),
    name = "-log(p)",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_side = "left",
    column_title = "Pain DEG Enrichment",
    heatmap_legend_param = list(
        title = "-log(p)",
        title_position = "topcenter", # Corrected positioning
        title_gp = gpar(fontsize = 10),
        labels_gp = gpar(fontsize = 8)
    ),
    column_names_side = "bottom",
    border = TRUE,  # Add a border around the heatmap
    border_gp = gpar(col = "black", lwd = 2)
)

# Display the heatmap
pdf(here("plots", "21_pain_enrichment", "single_nucleus_heatmap_up_down_adjpval_0.1.pdf"), heigh = 4, width = 2.5)
draw(heatmap_combined, merge_legend = F, annotation_legend_side = "bottom")
grid.text("",
          x = unit(0.5, "npc"), y = unit(0.02, "npc"),
          just = "center", gp = gpar(fontsize = 10, col = "black"))
dev.off()
