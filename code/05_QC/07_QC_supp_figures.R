setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("ggspavis"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("here"))

load(here("processed-data", "05_QC", "spe_QC.Rdata"), verbose = TRUE)
spe$discard_auto_br <- spe$low_sum_br | spe$low_detected_br
samples <- unique(spe$sample_id)
spe$discard_extreme <- spe$extreme_low_sum | spe$extreme_low_detected

# create 2-3 rows of vis_clus plots that show the "naive" filters for sum/det
p1 <- vis_clus(spe = spe, sampleid = samples[1], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    # remove legend because of formatting bug
    theme(legend.position="none")
p2 <- vis_clus(spe = spe, sampleid = samples[2], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p3 <- vis_clus(spe = spe, sampleid = samples[3], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1)  +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p4 <- vis_clus(spe = spe, sampleid = samples[4], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p5 <- vis_clus(spe = spe, sampleid = samples[5], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p6 <- vis_clus(spe = spe, sampleid = samples[6], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p7 <- vis_clus(spe = spe, sampleid = samples[7], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p8 <- vis_clus(spe = spe, sampleid = samples[8], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p9 <- vis_clus(spe = spe, sampleid = samples[9], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p10 <- vis_clus(spe = spe, sampleid = samples[10], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p11 <- vis_clus(spe = spe, sampleid = samples[11], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p12 <- vis_clus(spe = spe, sampleid = samples[12], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p13 <- vis_clus(spe = spe, sampleid = samples[13], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p14 <- vis_clus(spe = spe, sampleid = samples[14], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p15 <- vis_clus(spe = spe, sampleid = samples[15], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p16 <- vis_clus(spe = spe, sampleid = samples[16], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p17 <- vis_clus(spe = spe, sampleid = samples[17], clustervar = "discard_auto_br", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))

png(here("plots", "05_QC", "QC_discard_lowFeatures.png"), width = 25, height = 25, units = "in", res=300)

wrap_plots(p1, p2, p3, p4, p5, p6,
           p7, p8, p9, p10, p11, p12,
           p13, p14, p15, p16, p17,
           guides = "collect",
           nrow = 5) +
    plot_annotation(title="3 MADs Filter for Low Sum and Low Detected",theme=theme(plot.title = element_text(size = 40)))

dev.off()


# create 2-3 rows of vis_clus plots that show the "naive" filter for mito
p1 <- vis_gene(spe = spe, sampleid = samples[1], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[1])
p2 <- vis_gene(spe = spe, sampleid = samples[2], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[2])
p3 <- vis_gene(spe = spe, sampleid = samples[3], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1)  +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[3])
p4 <- vis_gene(spe = spe, sampleid = samples[4], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[4])
p5 <- vis_gene(spe = spe, sampleid = samples[5], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[5])
p6 <- vis_gene(spe = spe, sampleid = samples[6], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[6])
p7 <- vis_gene(spe = spe, sampleid = samples[7], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[7])
p8 <- vis_gene(spe = spe, sampleid = samples[8], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[8])
p9 <- vis_gene(spe = spe, sampleid = samples[9], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[9])
p10 <- vis_gene(spe = spe, sampleid = samples[10], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[10])
p11 <- vis_gene(spe = spe, sampleid = samples[11], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[11])
p12 <- vis_gene(spe = spe, sampleid = samples[12], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[12])
p13 <- vis_gene(spe = spe, sampleid = samples[13], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[13])
p14 <- vis_gene(spe = spe, sampleid = samples[14], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[14])
p15 <- vis_gene(spe = spe, sampleid = samples[15], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[15])
p16 <- vis_gene(spe = spe, sampleid = samples[16], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[16])
p17 <- vis_gene(spe = spe, sampleid = samples[17], assayname = "counts", geneid = "subsets_Mito_percent",  point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    ggtitle(samples[17])

png(here("plots", "05_QC", "QC_discard_highMito.png"), width = 25, height = 25, units = "in", res=300)

wrap_plots(p1, p2, p3, p4, p5, p6,
           p7, p8, p9, p10, p11, p12,
           p13, p14, p15, p16, p17,
           nrow = 5) +
    plot_annotation(title="Mitochondrial Percent",theme=theme(plot.title = element_text(size = 40)))

dev.off()

# create 2-3 rows of vis_clus plots that show the extreme filters we actually used
p1 <- vis_clus(spe = spe, sampleid = samples[1], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    # remove legend because of formatting bug
    theme(legend.position="none")
p2 <- vis_clus(spe = spe, sampleid = samples[2], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20)) +
    theme(legend.position="none")
p3 <- vis_clus(spe = spe, sampleid = samples[3], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1)  +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))  +
    theme(legend.position="none")
p4 <- vis_clus(spe = spe, sampleid = samples[4], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))  +
    theme(legend.position="none")
p5 <- vis_clus(spe = spe, sampleid = samples[5], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))  +
    theme(legend.position="none")
p6 <- vis_clus(spe = spe, sampleid = samples[6], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))  +
    theme(legend.position="none")
p7 <- vis_clus(spe = spe, sampleid = samples[7], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p8 <- vis_clus(spe = spe, sampleid = samples[8], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p9 <- vis_clus(spe = spe, sampleid = samples[9], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p10 <- vis_clus(spe = spe, sampleid = samples[10], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))  +
    theme(legend.position="none")
p11 <- vis_clus(spe = spe, sampleid = samples[11], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))  +
    theme(legend.position="none")
p12 <- vis_clus(spe = spe, sampleid = samples[12], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p13 <- vis_clus(spe = spe, sampleid = samples[13], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))  +
    theme(legend.position="none")
p14 <- vis_clus(spe = spe, sampleid = samples[14], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))  +
    theme(legend.position="none")
p15 <- vis_clus(spe = spe, sampleid = samples[15], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p16 <- vis_clus(spe = spe, sampleid = samples[16], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))
p17 <- vis_clus(spe = spe, sampleid = samples[17], clustervar = "discard_extreme", colors = c("FALSE" = "yellow", "TRUE" = "blue"), point_size = 1) +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(legend.text = element_text(size = 20))  +
    theme(legend.position="none")

png(here("plots", "05_QC", "QC_discard_extremeFeatures.png"), width = 25, height = 25, units = "in", res=300)

wrap_plots(p1, p2, p3, p4, p5, p6,
           p7, p8, p9, p10, p11, p12,
           p13, p14, p15, p16, p17,
           guides = "collect",
           nrow = 5) +
    plot_annotation(title="Extreme Filter for Sum < 20 and Detected < 20",theme=theme(plot.title = element_text(size = 40)))

dev.off()

# create violin plots of true/false for this new filter
