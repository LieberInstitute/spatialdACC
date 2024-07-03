# http://research.libd.org/DeconvoBuddies/articles/Deconvolution_with_Bisque.html

setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("spatialLIBD")
library("DeconvoBuddies")
library("SummarizedExperiment")
library("SingleCellExperiment")
library("BisqueRNA")
library("dplyr")
library("tidyr")
library("tibble")
library("here")
library("scran")

# Load the bulk data
load(here("processed-data", "PTSD_bulk", "rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata"))
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

# load the sce object
load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))
sce <- logNormCounts(sce)
rownames(sce) <- rowData(sce)$gene_id

marker_stats <- get_mean_ratio2(sce, cellType_col = "cellType_azimuth")

pdf(here("plots", "19_bulk_deconvolution", "marker_express_Astro.pdf"))
DeconvoBuddies::plot_marker_express(sce,
                                    stats = marker_stats,
                                    cell_type = "Astro",
                                    cellType_col = "cellType_azimuth",
                                    n_genes = 10,
                                    rank_col = "rank_ratio",
                                    anno_col = "anno_ratio",
)
dev.off()

marker_genes <- marker_stats |>
    filter(rank_ratio <= 25, gene %in% rownames(rse_gene)) |>
    pull(gene)

marker_genes <- unique(marker_genes)

length(marker_genes)
# 462

exp_set_bulk <- Biobase::ExpressionSet(
    assayData = assays(rse_gene[marker_genes, ])$counts,
    phenoData = AnnotatedDataFrame(
        as.data.frame(colData(rse_gene))[c("SAMPLE_ID")]
    )
)

exp_set_sce <- Biobase::ExpressionSet(
    assayData = as.matrix(assays(sce[marker_genes, ])$counts),
    phenoData = AnnotatedDataFrame(
        as.data.frame(colData(sce))[, c("cellType_azimuth", "brain")]
    )
)

## check for nuclei with 0 marker expression
zero_cell_filter <- colSums(exprs(exp_set_sce)) != 0
message("Exclude ", sum(!zero_cell_filter), " cells")
# Exclude 0 cells

exp_set_sce <- exp_set_sce[, zero_cell_filter]

est_prop <- ReferenceBasedDecomposition(
    bulk.eset = exp_set_bulk,
    sc.eset = exp_set_sce,
    cell.types = "cellType_azimuth",
    subject.names = "brain",
    use.overlap = FALSE
)

est_prop$bulk.props <- t(est_prop$bulk.props)

pd <- colData(rse_gene) |>
    as.data.frame() |>
    select(Sample = RNum, Sex, Group)

## make proportion estimates long so they are ggplot friendly
prop_long <- est_prop$bulk.props |>
    as.data.frame() |>
    tibble::rownames_to_column("Sample") |>
    tidyr::pivot_longer(!Sample, names_to = "cell_type", values_to = "prop") |>
    left_join(pd)

pdf(here("plots", "19_bulk_deconvolution", "single_nucleus_bulk_deconvolution_bisque.pdf"))
plot_composition_bar(prop_long = prop_long, sample_col = "Sample", x_col = "Group", min_prop_text = 0.025)
plot_composition_bar(prop_long = prop_long, sample_col = "Sample", x_col = "Sample", add_text = FALSE)
dev.off()
