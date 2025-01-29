setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library(here)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(EnhancedVolcano)

# Load DEG data
load(
    file = here("processed-data", "11_differential_expression", "pseudobulk", "nnSVG_precast_DE", "nnSVG_PRECAST_captureArea_9.Rdata")
)

# Subset modeling results
enrichment_results <- modeling_results[["enrichment"]]
dim(enrichment_results)
# 13576    30

# remove 1 duplicate
enrichment_results <- enrichment_results[!duplicated(enrichment_results$gene),]

# Load the bulk data
file = here("processed-data", "pain", "GSE227159_Processed-data_Becker-Lutz-Yalcin_BLA-ACC-pathway.csv")
bulk <- read.csv(file, row.names = 2)

dim(bulk)
# [1] 37786    50


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
bulk$p <- as.numeric(as.character(bulk$p))
bulk$adj_p <- as.numeric(as.character(bulk$adj_p))
bulk$logFC <- as.numeric(as.character(bulk$logFC))

orthology<-read.csv(file=here::here('raw-data','Retro-seq',
                                    'human_mouse_orthologs.csv'))

names <- orthology[orthology$Column3 %in% bulk$Gene.name,]

names <- names[match(bulk$Gene.name, names$Column3),]

setdiff(names$Column3, bulk$Gene.name)

sum(is.na(names$Column1))
# [1] 1436

bulk$Human.gene.name <- names$Column1

# 3 of the 54 DE genes do not have matching names

# remove non-matching human genes
bulk <- bulk[!is.na(bulk$Human.gene.name),]
dim(bulk)
# 12820 4

length(unique(bulk$Human.gene.name))
# [1] 12659

bulk[duplicated(bulk$Human.gene.name),]
# these genes are not of interest since they are not significant DEGs,
# just remove one instance of each to avoid errors downstream

bulk <- bulk[!duplicated(bulk$Human.gene.name),]

# Overlap to get genes in both datasets
overlap <- intersect(bulk$Human.gene.name, enrichment_results$gene)
bulk <- bulk[bulk$Human.gene.name %in% overlap, ]
enrichment_results <- enrichment_results[enrichment_results$gene %in% overlap, ]

sum(bulk$adj_p<0.05) # 144

    # make a volcano plot of the pain dataset to figure out a good threshold
    #volcano plots
    thresh_fdr <- 0.1
    thresh_logfc <- log2(1.25)

    fdrs <- bulk$adj_p
    logfc <- bulk$logFC

    sig <- (fdrs < thresh_fdr) & (abs(logfc) > thresh_logfc)
    print(table(sig))

    df_list <- data.frame(
        gene_name = bulk$Human.gene.name,
        logFC = logfc,
        FDR = fdrs,
        sig = sig
    )

    pdf(file = here::here("plots", "21_pain_enrichment", "pain_volcano.pdf"),
        width = 8.5, height = 8)

    print(EnhancedVolcano(df_list,
                          lab = df_list$gene_name,
                          x = 'logFC',
                          y = 'FDR',
                          FCcutoff = 0.25,
                          pCutoff = 0.1,
                          ylab = "-log10 FDR",
                          legendLabels = c('Not sig.','Log (base 2) FC','FDR',
                                           'FDR & Log (base 2) FC'),
                          title = "pain dataset",
                          subtitle = "after ortho matching and conversion to human gene names",
    )
    )

    dev.off()

generate_spatial_heatmap <- function(adj_pval_threshold = 0.1) {

    # Define upregulated and downregulated DEGs in bulk
    DE_bulk_up <- bulk[bulk$adj_p < adj_pval_threshold & bulk$logFC > 0, ]$Human.gene.name
    DE_bulk_down <- bulk[bulk$adj_p < adj_pval_threshold & bulk$logFC < 0, ]$Human.gene.name

    nonDE_bulk_up <- setdiff(bulk$Human.gene.name, DE_bulk_up)
    nonDE_bulk_down <- setdiff(bulk$Human.gene.name, DE_bulk_down)

    #p_val_threshold <- 0.05

    #DE_bulk_up <- bulk[bulk$p < p_val_threshold & bulk$logFC > 0, ]$Human.gene.name
    #DE_bulk_down <- bulk[bulk$p < p_val_threshold & bulk$logFC < 0, ]$Human.gene.name

    #nonDE_bulk_up <- setdiff(bulk$Human.gene.name, DE_bulk_up)
    #nonDE_bulk_down <- setdiff(bulk$Human.gene.name, DE_bulk_down)

    results_up <- list()
    results_down <- list()

    for (i in c("L1","L2","L3","L5","L6a","L6b","WM")) {
        #top_n <- 100

        #t_stat_threshold <- sort(enrichment_results[[paste0("t_stat_", i)]], decreasing = T)[top_n]

        #DE_clust_genes_up <- enrichment_results[enrichment_results[[paste0("t_stat_", i)]] >= t_stat_threshold, ]$gene
        #DE_clust_genes_down <- DE_clust_genes_up

        #nonDE_clust_genes_up <- enrichment_results[enrichment_results[[paste0("t_stat_", i)]] < t_stat_threshold, ]$gene
        #nonDE_clust_genes_down <- nonDE_clust_genes_up

        DE_clust_genes_up <- enrichment_results[
            enrichment_results[[paste0("fdr_", i)]] < 0.05 & enrichment_results[[paste0("logFC_", i)]] > 1, ]$gene
        DE_clust_genes_down <- DE_clust_genes_up

        print(i)
        print(length(DE_clust_genes_up))

        nonDE_clust_genes_up <- enrichment_results[
            enrichment_results[[paste0("fdr_", i)]] >= 0.05 | enrichment_results[[paste0("logFC_", i)]] <= 1, ]$gene
        nonDE_clust_genes_down <- nonDE_clust_genes_up

        DE_clust_DE_bulk_up <- length(intersect(DE_clust_genes_up, DE_bulk_up))
        DE_clust_nonDE_bulk_up <- length(intersect(DE_clust_genes_up, nonDE_bulk_up))
        nonDE_clust_DE_bulk_up <- length(intersect(nonDE_clust_genes_up, DE_bulk_up))
        nonDE_clust_nonDE_bulk_up <- length(intersect(nonDE_clust_genes_up, nonDE_bulk_up))

        contingency_table_up <- matrix(c(DE_clust_DE_bulk_up, DE_clust_nonDE_bulk_up,
                                         nonDE_clust_DE_bulk_up, nonDE_clust_nonDE_bulk_up),
                                       nrow = 2, byrow = TRUE,
                                       dimnames = list("Cluster" = c("DE", "nonDE"),
                                                       "Bulk" = c("DE", "nonDE")))

        fisher_test_up <- fisher.test(contingency_table_up)
        results_up[[i]] <- list(fisher_test = fisher_test_up)

        DE_clust_DE_bulk_down <- length(intersect(DE_clust_genes_down, DE_bulk_down))
        DE_clust_nonDE_bulk_down <- length(intersect(DE_clust_genes_down, nonDE_bulk_down))
        nonDE_clust_DE_bulk_down <- length(intersect(nonDE_clust_genes_down, DE_bulk_down))
        nonDE_clust_nonDE_bulk_down <- length(intersect(nonDE_clust_genes_down, nonDE_bulk_down))

        contingency_table_down <- matrix(c(DE_clust_DE_bulk_down, DE_clust_nonDE_bulk_down,
                                           nonDE_clust_DE_bulk_down, nonDE_clust_nonDE_bulk_down),
                                         nrow = 2, byrow = TRUE,
                                         dimnames = list("Cluster" = c("DE", "nonDE"),
                                                         "Bulk" = c("DE", "nonDE")))

        fisher_test_down <- fisher.test(contingency_table_down)
        results_down[[i]] <- list(fisher_test = fisher_test_down)
    }

    pvalues_up <- sapply(results_up, function(x) x$fisher_test$p.value)
    pvalues_down <- sapply(results_down, function(x) x$fisher_test$p.value)

    combined_pvalues <- cbind(
        Upreg = pvalues_up,
        Downreg = pvalues_down
    )

    row_means <- rowMeans(-log10(combined_pvalues), na.rm = TRUE)
    col_means <- colMeans(-log10(combined_pvalues), na.rm = TRUE)

    ordered_rows <- order(row_means, decreasing = TRUE)
    ordered_cols <- order(col_means, decreasing = TRUE)

    combined_pvalues_ordered <- combined_pvalues[ordered_rows, ordered_cols]

    combined_pvalues_ordered <- combined_pvalues_ordered[c(1,3,2,4,5,6,7),]

    col_fun <- colorRamp2(
        c(1.3, max(-log10(combined_pvalues_ordered))),
        c("white", "blue") # White for -log10(p) >= 1.3 (p >= 0.05), blue for more significant p-values
    )

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
            title_position = "topcenter",
            title_gp = gpar(fontsize = 10),
            labels_gp = gpar(fontsize = 8)
        ),
        column_names_side = "bottom",
        border = TRUE,
        border_gp = gpar(col = "black", lwd = 2)
    )

    print(combined_pvalues_ordered)

    # Return the heatmap object
    return(heatmap_combined)
}

# Display the heatmap
pdf(here("plots", "21_pain_enrichment", "spatial_heatmap_up_down_adjpval_0.1.pdf"),
    height = 4, width = 2)
generate_spatial_heatmap()
grid.text("",
          x = unit(0.5, "npc"), y = unit(0.02, "npc"),
          just = "center", gp = gpar(fontsize = 10, col = "black"))
dev.off()
