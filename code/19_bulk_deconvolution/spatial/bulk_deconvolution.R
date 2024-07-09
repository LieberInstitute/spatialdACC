setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library("spatialLIBD")
library("DeconvoBuddies")
library("SummarizedExperiment")
library("SpatialExperiment")
library("BisqueRNA")
library("dplyr")
library("tidyr")
library("tibble")
library("here")
library("scran")

# Load the bulk data
load(here("processed-data", "PTSD_bulk", "rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata"))
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

# load the spe object
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9.Rdata"))
spe$PRECAST_cluster <- unfactor(spe$PRECAST_cluster)
spe$PRECAST_cluster[spe$PRECAST_cluster == 3] <- "WM1"
spe$PRECAST_cluster[spe$PRECAST_cluster == 8] <- "WM2"
spe$PRECAST_cluster[spe$PRECAST_cluster == 7] <- "WM-CC"
spe$PRECAST_cluster[spe$PRECAST_cluster == 5] <- "L6b"
spe$PRECAST_cluster[spe$PRECAST_cluster == 6] <- "L6a"
spe$PRECAST_cluster[spe$PRECAST_cluster == 4] <- "L5"
spe$PRECAST_cluster[spe$PRECAST_cluster == 2] <- "L3"
spe$PRECAST_cluster[spe$PRECAST_cluster == 1] <- "L2"
spe$PRECAST_cluster[spe$PRECAST_cluster == 9] <- "L1"

marker_stats <- get_mean_ratio2(spe, cellType_col = "PRECAST_cluster")

marker_genes <- marker_stats |>
    filter(rank_ratio <= 25, gene %in% rownames(rse_gene)) |>
    pull(gene)

marker_genes <- unique(marker_genes)

length(marker_genes)
# 187

exp_set_bulk <- Biobase::ExpressionSet(
    assayData = assays(rse_gene[marker_genes, ])$counts,
    phenoData = AnnotatedDataFrame(
        as.data.frame(colData(rse_gene))[c("SAMPLE_ID")]
    )
)

exp_set_spe <- Biobase::ExpressionSet(
    assayData = as.matrix(assays(spe[marker_genes, ])$counts),
    phenoData = AnnotatedDataFrame(
        as.data.frame(colData(spe))[, c("PRECAST_cluster", "brnum")]
    )
)

## check for nuclei with 0 marker expression
zero_cell_filter <- colSums(exprs(exp_set_spe)) != 0
message("Exclude ", sum(!zero_cell_filter), " cells")
# Exclude 0 cells

exp_set_spe <- exp_set_spe[, zero_cell_filter]

est_prop <- ReferenceBasedDecomposition(
    bulk.eset = exp_set_bulk,
    sc.eset = exp_set_spe,
    cell.types = "PRECAST_cluster",
    subject.names = "brnum",
    use.overlap = FALSE
)

est_prop$bulk.props <- t(est_prop$bulk.props)

pd <- colData(rse_gene) |>
    as.data.frame() |>
    select(Sample = RNum, Sex, Group)

# create a new label for Group
# make PTSD and MDD the same group called "PTSD or MDD"
pd$Group <- ifelse(pd$Group %in% c("PTSD", "MDD"), "PTSD or MDD", pd$Group)

## make proportion estimates long so they are ggplot friendly
prop_long <- est_prop$bulk.props |>
    as.data.frame() |>
    tibble::rownames_to_column("Sample") |>
    tidyr::pivot_longer(!Sample, names_to = "cell_type", values_to = "prop") |>
    left_join(pd)

pdf(here("plots", "19_bulk_deconvolution", "spatial_bulk_deconvolution_bisque.pdf"))
plot_composition_bar(prop_long = prop_long, sample_col = "Sample", x_col = "Group", min_prop_text = 0.025)
plot_composition_bar(prop_long = prop_long, sample_col = "Sample", x_col = "Sample", add_text = FALSE)
dev.off()
