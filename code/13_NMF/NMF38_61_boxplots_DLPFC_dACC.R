setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(SpatialExperiment)
library(scater)
library(RcppML)
library(ggspavis)
library(here)
library(scRNAseq)
library(Matrix)
library(scran)
library(scuttle)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(igraph)
library(bluster)
library(patchwork)
library(cowplot)
library(projectR)
library(spatialLIBD)
library(gridExtra)
library(escheR)
library(dplyr)
library(tidyverse)

load(file = here("processed-data", "13_NMF", "spe_NMF.Rdata"))
spe_dACC <- spe

load(file = here("processed-data", "13_NMF", "spe_NMF_DLPFC_30.Rdata"))
spe_DLPFC <- spe

dat_dACC <- as.data.frame(colData(spe_dACC))
dat_DLPFC <- as.data.frame(colData(spe_DLPFC))

# only keep L5
dat_dACC <- dat_dACC[which(dat_dACC$layer == "L5"),]
summary_dACC <- dat_dACC %>%
    group_by(brnum) %>%
    summarize(totalL5_dACC = n(), count38_dACC = sum(nmf38 > 0), count61_dACC = sum(nmf61 > 0))

summary_dACC$frac38_dACC <- summary_dACC$count38_dACC / summary_dACC$totalL5_dACC
summary_dACC$frac61_dACC <- summary_dACC$count61_dACC / summary_dACC$totalL5_dACC

summary_dACC$region <- rep("dACC", 10)

# only keep L5
dat_DLPFC <- dat_DLPFC[which(dat_DLPFC$BayesSpace_harmony_09 == 4),]
summary_DLPFC <- dat_DLPFC %>%
    group_by(subject) %>%
    summarize(totalL5_DLPFC = n(), count38_DLPFC = sum(nmf38 > 0), count61_DLPFC = sum(nmf61 > 0))

summary_DLPFC$frac38_DLPFC <- summary_DLPFC$count38_DLPFC / summary_DLPFC$totalL5_DLPFC
summary_DLPFC$frac61_DLPFC <- summary_DLPFC$count61_DLPFC / summary_DLPFC$totalL5_DLPFC

summary_DLPFC$region <- rep("dlPFC", 10)

summary_overall <- summary_dACC
summary_overall[c(11:20),] <- summary_DLPFC

p1 <- ggplot(summary_overall, aes(x=region, y=frac38_dACC)) +
    geom_boxplot(outlier.shape = NA, color="#FFD700") +
    ylim(c(0,0.4)) +
    geom_point(color="#FFD700", size=1, alpha=0.8) +
    ylab("Frac. Nonzero NMF38") +
    ggtitle("") +
    theme_bw()

p2 <- ggplot(summary_overall, aes(x=region, y=frac61_dACC)) +
    geom_boxplot(outlier.shape = NA, color="#FFD700") +
    geom_point(color="#FFD700", size=1, alpha=0.8) +
    ylab("Frac. Nonzero NMF61") +
    ylim(c(0,0.4)) +
    ggtitle("") +
    theme_bw()

# also calculate average weights for SRT data
dat_dACC <- as.data.frame(colData(spe_dACC))
dat_DLPFC <- as.data.frame(colData(spe_DLPFC))

dat_DLPFC <- dat_DLPFC[which(dat_DLPFC$BayesSpace_harmony_09 == 4),]
summary_DLPFC <- dat_DLPFC %>%
    group_by(subject) %>%
    summarize(avg38 = mean(nmf38), avg61 = mean(nmf61))

summary_DLPFC$region <- rep("dlPFC", 10)

dat_dACC <- dat_dACC[which(dat_dACC$layer == "L5"),]
summary_dACC <- dat_dACC %>%
    group_by(brnum) %>%
    summarize(avg38 = mean(nmf38), avg61 = mean(nmf61))

summary_dACC$region <- rep("dACC", 10)

summary_overall <- summary_dACC
summary_overall[c(11:20),] <- summary_DLPFC


p3 <- ggplot(summary_overall, aes(x=region, y=avg38, color=region)) +
    geom_boxplot(outlier.shape = NA) +
    ylim(c(0,0.00011)) +
    geom_point(size=1, alpha=0.8) +
    ylab("Avg. NMF38 Weight") +
    ggtitle("") +
    xlab("region") +
    theme_bw() +
    scale_color_manual(values = c("#FFD700", "#FFD700")) +
    theme(legend.position="none")

p4 <- ggplot(summary_overall, aes(x=region, y=avg61, color=region)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(size=1, alpha=0.8) +
    ylab("Avg. NMF61 Weight") +
    ylim(c(0.000,0.000185)) +
    ggtitle("") +
    xlab("region") +
    theme_bw() +
    scale_color_manual(values = c("#FFD700", "#FFD700")) +
    theme(legend.position="none")


# calculate nonzero spots for snRNA-seq data
load(file = here("processed-data", "13_NMF", "DLPFC_dACC_celltype_NMF.Rdata"))

dat_DLPFC_orig <- dat_DLPFC

dat_DLPFC <- dat_DLPFC[which(dat_DLPFC$cellType_layer %in% c("Excit_L5")),]
summary_DLPFC <- dat_DLPFC %>%
    group_by(BrNum) %>%
    summarize(totalL5_DLPFC = n(), count38_DLPFC = sum(nmf38 > 0), count61_DLPFC = sum(nmf61 > 0))

summary_DLPFC$region <- rep("dlPFC", dim(summary_DLPFC)[1])
summary_DLPFC$celltype <- rep("Excit_L5", dim(summary_DLPFC)[1])

summary_DLPFC$frac38_DLPFC <- summary_DLPFC$count38_DLPFC / summary_DLPFC$totalL5_DLPFC
summary_DLPFC$frac61_DLPFC <- summary_DLPFC$count61_DLPFC / summary_DLPFC$totalL5_DLPFC

dat_dACC_orig <- dat_dACC
dat_dACC <- dat_dACC[which(dat_dACC$cellType_azimuth %in% c("L5_IT")),]
summary_dACC_IT <- dat_dACC %>%
    group_by(brain) %>%
    summarize(totalL5_dACC = n(), count38_dACC = sum(NMF_38 > 0), count61_dACC = sum(NMF_61 > 0))

summary_dACC_IT$region <- rep("dACC", 10)
summary_dACC_IT$celltype <- rep("L5_IT", 10)

summary_dACC_IT$frac38_dACC <- summary_dACC_IT$count38_dACC / summary_dACC_IT$totalL5_dACC
summary_dACC_IT$frac61_dACC <- summary_dACC_IT$count61_dACC / summary_dACC_IT$totalL5_dACC


dat_dACC <- dat_dACC_orig
dat_dACC <- dat_dACC[which(dat_dACC$cellType_azimuth %in% c("L5_ET")),]
summary_dACC_ET <- dat_dACC %>%
    group_by(brain) %>%
    summarize(totalL5_dACC = n(), count38_dACC = sum(NMF_38 > 0), count61_dACC = sum(NMF_61 > 0))

summary_dACC_ET$region <- rep("dACC", 10)
summary_dACC_ET$celltype <- rep("L5_ET", 10)

summary_dACC_ET$frac38_dACC <- summary_dACC_ET$count38_dACC / summary_dACC_ET$totalL5_dACC
summary_dACC_ET$frac61_dACC <- summary_dACC_ET$count61_dACC / summary_dACC_ET$totalL5_dACC

summary_overall <- summary_dACC_IT
summary_overall[c(11:20),] <- summary_DLPFC

p5 <- ggplot(summary_overall, aes(x=region, y=frac38_dACC, color=region)) +
    geom_boxplot(outlier.shape = NA) +
    ylim(c(0,1)) +
    geom_point(size=1, alpha=0.8) +
    ylab("Frac. Nonzero NMF38") +
    scale_color_manual(values = c("#6A33C2", "black")) +
    ggtitle("") +
    theme_bw() +
    theme(legend.position="none")

summary_overall <- summary_dACC_ET
summary_overall[c(11:20),] <- summary_DLPFC

p6 <- ggplot(summary_overall, aes(x=region, y=frac61_dACC, color=region)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(size=1, alpha=0.8) +
    ylab("Frac. Nonzero NMF61") +
    ylim(c(0,1)) +
    ggtitle("") +
    theme_bw() +
    scale_color_manual(values = c("#C8308C", "black")) +
    theme(legend.position="none")


# calculate average weights for snRNA-seq data
load(file = here("processed-data", "13_NMF", "DLPFC_dACC_celltype_NMF.Rdata"))

dat_DLPFC_orig <- dat_DLPFC

dat_DLPFC <- dat_DLPFC[which(dat_DLPFC$cellType_layer %in% c("Excit_L5")),]
summary_DLPFC <- dat_DLPFC %>%
    group_by(BrNum) %>%
    summarize(avg38 = mean(nmf38), avg61 = mean(nmf61))

summary_DLPFC$region <- rep("dlPFC", dim(summary_DLPFC)[1])
summary_DLPFC$celltype <- rep("Excit_L5", dim(summary_DLPFC)[1])

dat_dACC_orig <- dat_dACC
dat_dACC <- dat_dACC[which(dat_dACC$cellType_azimuth %in% c("L5_IT")),]
summary_dACC_IT <- dat_dACC %>%
    group_by(brain) %>%
    summarize(avg38 = mean(NMF_38), avg61 = mean(NMF_61))

summary_dACC_IT$region <- rep("dACC", 10)
summary_dACC_IT$celltype <- rep("L5_IT", 10)

dat_dACC <- dat_dACC_orig
dat_dACC <- dat_dACC[which(dat_dACC$cellType_azimuth %in% c("L5_ET")),]
summary_dACC_ET <- dat_dACC %>%
    group_by(brain) %>%
    summarize(avg38 = mean(NMF_38), avg61 = mean(NMF_61))

summary_dACC_ET$region <- rep("dACC", 10)
summary_dACC_ET$celltype <- rep("L5_ET", 10)

summary_overall <- summary_dACC_IT
summary_overall[c(11:20),] <- summary_DLPFC

p7 <- ggplot(summary_overall, aes(x=region, y=avg38, color=region)) +
    geom_boxplot(outlier.shape = NA) +
    ylim(c(0,0.00065)) +
    geom_point(size=1, alpha=0.8) +
    ylab("Avg. NMF38 Weight") +
    ggtitle("") +
    xlab("region") +
    theme_bw() +
    scale_color_manual(values = c("#6A33C2", "black")) +
    theme(legend.position="none")

summary_overall <- summary_dACC_ET
summary_overall[c(11:20),] <- summary_DLPFC

p8 <- ggplot(summary_overall, aes(x=region, y=avg61, color=region)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(size=1, alpha=0.8) +
    ylab("Avg. NMF61 Weight") +
    ylim(c(0.000,0.003)) +
    ggtitle("") +
    xlab("region") +
    theme_bw() +
    scale_color_manual(values = c("#C8308C", "black")) +
    theme(legend.position="none")


pdf(file = here::here("plots", "13_NMF", "NMF_boxplots_DLPFC_dACC.pdf"), height = 6, width = 6)
wrap_plots(list(p1,p3,p2,p4), nrow = 2) + plot_annotation(title = "NMF38 and 61 in SRT",
                                                          caption = "Layer 5 in dACC and dlPFC")
wrap_plots(list(p5,p7,p6,p8), nrow = 2) + plot_annotation(title = "NMF38 and 61 in snRNA-seq",
                                                          caption = "L5 in dlPFC; L5 IT in dACC (top), L5 ET in dACC (bottom)")
dev.off()

