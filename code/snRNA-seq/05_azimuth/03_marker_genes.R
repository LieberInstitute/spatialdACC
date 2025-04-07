setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SingleCellExperiment)
library(here)
library(scran)
library(ggplot2)
library(spatialLIBD)
library(pheatmap)
library(scater)
library(patchwork)

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))
sce <- logNormCounts(sce)

marker.info <- scoreMarkers(sce, colData(sce)$cellType_azimuth)
marker.info

save(marker.info, file = here("processed-data", "snRNA-seq", "05_azimuth", "marker_info_azimuth.Rdata"))

# get top 20 markers for each cluster
top20_markers <- lapply(marker.info, function(x) {
  x <- x[order(x$mean.AUC, decreasing = TRUE),]
  x <- x[1:20,]
  return(x)
})

# write to csv file
for (i in 1:length(top20_markers)) {
  write.csv(top20_markers[[i]], file = here("processed-data", "snRNA-seq", "05_azimuth", paste0("top20_markers_", names(top20_markers)[i], ".csv")))
}

# visualize some of the layer markers recovered in the snRNA-seq data
sce$counts_PCP4 <- counts(sce)[which(rowData(sce)$gene_name=="PCP4"),]
sce$logcounts_PCP4 <- logcounts(sce)[which(rowData(sce)$gene_name=="PCP4"),]

sce$counts_RORB <- counts(sce)[which(rowData(sce)$gene_name=="RORB"),]
sce$logcounts_RORB <- logcounts(sce)[which(rowData(sce)$gene_name=="RORB"),]

sce$counts_ADCYAP1 <- counts(sce)[which(rowData(sce)$gene_name=="ADCYAP1"),]
sce$logcounts_ADCYAP1 <- logcounts(sce)[which(rowData(sce)$gene_name=="ADCYAP1"),]

sce$counts_ARHGAP4 <- counts(sce)[which(rowData(sce)$gene_name=="ARHGAP4"),]
sce$logcounts_ARHGAP4 <- logcounts(sce)[which(rowData(sce)$gene_name=="ARHGAP4"),]

sce$counts_CPLX3 <- counts(sce)[which(rowData(sce)$gene_name=="CPLX3"),]
sce$logcounts_CPLX3 <- logcounts(sce)[which(rowData(sce)$gene_name=="CPLX3"),]

sce$counts_KCTD8 <- counts(sce)[which(rowData(sce)$gene_name=="KCTD8"),]
sce$logcounts_KCTD8 <- logcounts(sce)[which(rowData(sce)$gene_name=="KCTD8"),]

sce$counts_NXPH3 <- counts(sce)[which(rowData(sce)$gene_name=="NXPH3"),]
sce$logcounts_NXPH3 <- logcounts(sce)[which(rowData(sce)$gene_name=="NXPH3"),]

sce$counts_DRD5 <- counts(sce)[which(rowData(sce)$gene_name=="DRD5"),]
sce$logcounts_DRD5 <- logcounts(sce)[which(rowData(sce)$gene_name=="DRD5"),]

sce$counts_HBA1 <- counts(sce)[which(rowData(sce)$gene_name=="HBA1"),]
sce$logcounts_HBA1 <- logcounts(sce)[which(rowData(sce)$gene_name=="HBA1"),]

sce$counts_RSPO2 <- counts(sce)[which(rowData(sce)$gene_name=="RSPO2"),]
sce$logcounts_RSPO2 <- logcounts(sce)[which(rowData(sce)$gene_name=="RSPO2"),]

sce$counts_CD24 <- counts(sce)[which(rowData(sce)$gene_name=="CD24"),]
sce$logcounts_CD24 <- logcounts(sce)[which(rowData(sce)$gene_name=="CD24"),]

sce$counts_PVALB <- counts(sce)[which(rowData(sce)$gene_name=="PVALB"),]
sce$logcounts_PVALB <- logcounts(sce)[which(rowData(sce)$gene_name=="PVALB"),]

pdf(file = here::here("plots", "snRNA-seq", "05_azimuth", "marker_violin_plots.pdf"),
    width = 15, height = 6)

plotColData(sce, x = "cellType_azimuth", y = "counts_PCP4") +
    ggtitle("PCP4 Counts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "logcounts_PCP4") +
    ggtitle("PCP4 Logcounts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "counts_RORB") +
    ggtitle("RORB Counts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "logcounts_RORB") +
    ggtitle("RORB Logcounts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "counts_ADCYAP1") +
    ggtitle("ADCYAP1 Counts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "logcounts_ADCYAP1") +
    ggtitle("ADCYAP1 Logcounts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "counts_ARHGAP4") +
    ggtitle("ARHGAP4 Counts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "logcounts_ARHGAP4") +
    ggtitle("ARHGAP4 Logcounts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "counts_CPLX3") +
    ggtitle("CPLX3 Counts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "logcounts_CPLX3") +
    ggtitle("CPLX3 Logcounts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)


plotColData(sce, x = "cellType_azimuth", y = "counts_KCTD8") +
    ggtitle("KCTD8 Counts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "logcounts_KCTD8") +
    ggtitle("KCTD8 Logcounts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "counts_NXPH3") +
    ggtitle("NXPH3 Counts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "logcounts_NXPH3") +
    ggtitle("NXPH3 Logcounts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "counts_DRD5") +
    ggtitle("DRD5 Counts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "logcounts_DRD5") +
    ggtitle("DRD5 Logcounts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "counts_HBA1") +
    ggtitle("HBA1 Counts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "logcounts_HBA1") +
    ggtitle("HBA1 Logcounts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "counts_RSPO2") +
    ggtitle("RSPO2 Counts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "logcounts_RSPO2") +
    ggtitle("RSPO2 Logcounts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "counts_CD24") +
    ggtitle("CD24 Counts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "logcounts_CD24") +
    ggtitle("CD24 Logcounts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "counts_PVALB") +
    ggtitle("PVALB Counts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

plotColData(sce, x = "cellType_azimuth", y = "logcounts_PVALB") +
    ggtitle("PVALB Logcounts by Azimuth Cluster") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1)

dev.off()



# visualize 3 of the layer markers recovered in the snRNA-seq data
load(file = here("processed-data", "snRNA-seq", "05_azimuth", "celltype_colors.Rdata"))

#remove Sst chodl
idx <- which(sce$cellType_azimuth=="Sst Chodl")
sce <- sce[-idx,]
sce <- sce[,-idx]

sce$logcounts_PCP4 <- logcounts(sce)[which(rowData(sce)$gene_name=="PCP4"),]

sce$logcounts_CPLX3 <- logcounts(sce)[which(rowData(sce)$gene_name=="CPLX3"),]

sce$logcounts_CD24 <- logcounts(sce)[which(rowData(sce)$gene_name=="CD24"),]

p1 <- plotColData(sce, x = "cellType_azimuth", y = "logcounts_PCP4", colour_by = "cellType_azimuth") +
    ggtitle("Marker Gene Expression by Azimuth Cell Type") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1) +
    scale_color_manual(values = celltype_colors) +
    theme(legend.position="none") +
    xlab("") +
    ylab("PCP4 Logcounts") +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
    )

p2 <- plotColData(sce, x = "cellType_azimuth", y = "logcounts_CPLX3", colour_by = "cellType_azimuth") +
    ggtitle("") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1) +
    scale_color_manual(values = celltype_colors) +
    theme(legend.position="none") +
    xlab("") +
    ylab("CPLX3 Logcounts") +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
    )

p3 <- plotColData(sce, x = "cellType_azimuth", y = "logcounts_CD24", colour_by = "cellType_azimuth") +
    ggtitle("") +
    facet_wrap(~ sce$cellType_azimuth, scales = "free_x", nrow = 1) +
    scale_color_manual(values = celltype_colors) +
    theme(legend.position="none") +
    xlab("") +
    ylab("CD24 Logcounts") +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
    )


pdf(file = here::here("plots", "snRNA-seq", "05_azimuth", "subset_marker_violin_plots.pdf"),
    width = 15, height = 7)

wrap_plots(p1,p2,p3, nrow=3)

dev.off()

