setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("here")
library("spatialLIBD")
library("SpatialExperiment")
library("tidyverse")
library("ggplot2")
library("scran")

# plot pseuobulked samples for WM and CC in a scatterplot maybe
load(file = here("processed-data", "20_WM_comparisons", "pseudobulk_manual_anno.RData"))
logcounts_matrix <- logcounts(spe_pseudo)

sample_names <- colnames(logcounts_matrix)
donor <- substr(sample_names, 1, 13)
region <- sapply(strsplit(sample_names, "_"), function(x) tail(x, 1))

genes_of_interest <- c("TTR", "MSX1", "MSX2", "NTS")
regions_of_interest <- c("WM", "CC")

logcounts_matrix_filtered <- logcounts_matrix[rowData(spe_pseudo)$gene_name %in% genes_of_interest, ]
logcounts_matrix_filtered <- logcounts_matrix_filtered[, region %in% regions_of_interest]

df_plot <- data.frame(
    logcounts = as.vector(logcounts_matrix_filtered),
    gene = rep(rownames(logcounts_matrix_filtered), times = ncol(logcounts_matrix_filtered)),
    donor = rep(donor[region %in% regions_of_interest], each = nrow(logcounts_matrix_filtered)),
    region = rep(region[region %in% regions_of_interest], each = nrow(logcounts_matrix_filtered))
)

# replace gene id with gene name
df_plot$gene <- rowData(spe_pseudo)$gene_name[match(df_plot$gene, rowData(spe_pseudo)$gene_id)]

p <- ggplot(df_plot, aes(x = gene, y = logcounts, color = donor, shape = region)) +
    geom_point(position = position_jitter(width = 0.2, height = 0)) +  # Add jitter for better visibility
    theme_minimal() +
    labs(title = "Pseudobulked Logcounts in WM and CC",
         x = "Gene",
         y = "Logcounts") +
    scale_color_manual(values = RColorBrewer::brewer.pal(n = length(unique(df_plot$donor)), name = "Set2")) +  # Use color palette for donors
    scale_shape_manual(values = c("WM" = 16, "CC" = 17)) +  # Different shapes for regions
    theme(
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

# Save the plot as a PDF
pdf(file = here("plots", "20_WM_comparisons", "CP_WM_CC_genes_by_donor.pdf"), width = 6, height = 6)
print(p)
dev.off()
