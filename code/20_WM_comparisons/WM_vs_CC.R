setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("here")
library("spatialLIBD")
library("SpatialExperiment")
library("tidyverse")
library("escheR")
library("gridExtra")
library('EnhancedVolcano')
library("scater")

load(file = here("processed-data", "20_WM_comparisons", "spe_anno.Rdata"))
spe <- spe_anno

# vis MBP using escheR
spe$counts_MBP <- counts(spe)[which(rowData(spe)$gene_name=="MBP"),]

brains <- unique(spe$brnum)

pdf(file = here::here("plots", "20_WM_comparisons", "MBP_spot_plots_manual_anno.pdf"),
    width = 21, height = 20)

for (j in seq_along(brains)){
    speb <- spe[, which(spe$brnum == brains[j])]
    samples <- unique(speb$sample_id)
    print(length(samples))

    if (length(samples) == 1){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="counts_MBP", point_size = 9) |> add_ground(var="anno_label", stroke=0.5, point_size = 9) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        grid.arrange(p1, nrow = 1)
    } else if (length(samples) == 2){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="counts_MBP", point_size = 4) |> add_ground(var="anno_label", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="counts_MBP", point_size = 4) |> add_ground(var="anno_label", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="counts_MBP", point_size = 4) |> add_ground(var="anno_label", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="counts_MBP", point_size = 4) |> add_ground(var="anno_label", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))

        spe_3 <- speb[, which(speb$sample_id == samples[3])]
        p3 <- make_escheR(spe_3) |> add_fill(var="counts_MBP", point_size = 4) |> add_ground(var="anno_label", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[3]))
        grid.arrange(p1, p2, p3, nrow = 2)
    }
}

dev.off()

# exclude some samples from this comparison because the CC was not well detected
# V12Y31−080_B1, V12N28−334_C1, V12N28−334_B1
spe_anno <- spe_anno[, !spe_anno$sample_id %in% c("V12Y31-080_B1", "V12N28-334_C1", "V12N28-334_B1")]

# exclude one sample without CC
# V12N28−334_A1
spe_anno <- spe_anno[,!spe_anno$sample_id %in% c("V12N28-334_A1")]

dim(spe_anno)
#[1] 36601 27146

# relabel some clusters
# L2/3 = L2_3
# L4/5 = L4_5

spe_anno$anno_label <- gsub("L2/3", "L2_3", spe_anno$anno_label)
spe_anno$anno_label <- gsub("L4/5", "L4_5", spe_anno$anno_label)


spe_pseudo <-
    registration_pseudobulk(spe_anno,
                            var_registration = "anno_label",
                            var_sample_id = "sample_id",
                            min_ncells = 10
    )

dim(spe_pseudo)
#  15474    36

pca <- prcomp(t(assays(spe_pseudo)$logcounts))
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
pca_pseudo <- pca$x[, seq_len(20)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

## save pseudobulked spe file
save(
    spe_pseudo,
    file = here("processed-data", "20_WM_comparisons", "pseudobulk_manual_anno.RData")
)


## Plot PCs

pdf(file = here("plots", "20_WM_comparisons", "pseudobulk_anno.pdf"),
    width = 10, height = 10)

plotPCA(
    spe_pseudo,
    colour_by = "anno_label",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(spe_pseudo)$PCA_var_explained
)

plotPCA(
    spe_pseudo,
    colour_by = "sample_id",
    ncomponents = 2,
    point_size = 2,
    percentVar = metadata(spe_pseudo)$PCA_var_explained
)


vars <- getVarianceExplained(spe_pseudo,
                             variables = c("anno_label","sample_id")
)


plotExplanatoryVariables(vars)

dev.off()

modeling_results <- registration_wrapper(
    spe_anno,
    var_registration = "anno_label",
    var_sample_id = "sample_id",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

# save the modeling results
saveRDS(modeling_results, file = here("processed-data", "20_WM_comparisons", "modeling_results_WM_vs_CC.RDS"))

# find genes that are differentially expressed between WM and CC
# CC-WM such that CC is greater than WM
sig_genes <- sig_genes_extract(
    n = 20,
    modeling_results = modeling_results,
    model_type = "pairwise",
    reverse = FALSE
)[c(121:140),]

# save the significant genes
write.csv(sig_genes, file = here("processed-data", "20_WM_comparisons", "sig_genes_CC-WM.csv"))

# make volcano plots
thresh_fdr <- 0.05
thresh_logfc <- log2(1.5)
fdrs_gene_ids <- rowData(spe_anno)$gene_id
fdrs_gene_names <- rowData(spe_anno)$gene_name
fdrs <- modeling_results[["pairwise"]][,"fdr_CC-WM"]
logfc <- modeling_results[["pairwise"]][,"logFC_CC-WM"]

# Identify significant genes (low FDR and high logFC)
sig <- (fdrs < thresh_fdr) & (abs(logfc) > thresh_logfc)
# 271 genes

df_list <- data.frame(
    gene_name = modeling_results[["pairwise"]]$gene,
    logFC = logfc,
    FDR = fdrs,
    sig = sig
)

pdf(here("plots", "20_WM_comparisons", "volcano_plot_CC-WM.pdf"), width = 10, height = 10)
print(EnhancedVolcano(df_list,
                      lab = df_list$gene_name,
                      x = 'logFC',
                      y = 'FDR',
                      FCcutoff = 1.5,
                      pCutoff = 0.05,
                      ylab = "-log10 FDR",
                      legendLabels = c('Not sig.','Log (base 2) FC','FDR',
                                       'FDR & Log (base 2) FC'),
                      title = "manual annotations",
                      subtitle = paste0("CC vs. WM")
))
dev.off()
