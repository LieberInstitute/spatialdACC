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

load(file = here("processed-data", "13_NMF", "DLPFC_dACC_celltype_NMF.Rdata"))

dat_DLPFC <- as.data.frame(colData(sce_DLPFC))
dat_DLPFC <- dat_DLPFC[which(sce_DLPFC$cellType_azimuth %in% c("L5_IT")),]
summary_DLPFC_IT <- dat_DLPFC %>%
    group_by(Sample) %>%
    summarize(avg38 = mean(nmf38), avg61 = mean(nmf61))

summary_DLPFC_IT$region <- rep("dlPFC", dim(summary_DLPFC_IT)[1])
summary_DLPFC_IT$celltype <- rep("L5_IT", dim(summary_DLPFC_IT)[1])

dat_DLPFC <- as.data.frame(colData(sce_DLPFC))
dat_DLPFC <- dat_DLPFC[which(sce_DLPFC$cellType_azimuth %in% c("L5_ET")),]
summary_DLPFC_ET <- dat_DLPFC %>%
    group_by(Sample) %>%
    summarize(avg38 = mean(nmf38), avg61 = mean(nmf61))

summary_DLPFC_ET$region <- rep("dlPFC", dim(summary_DLPFC_ET)[1])
summary_DLPFC_ET$celltype <- rep("L5_ET", dim(summary_DLPFC_ET)[1])

dat_dACC <- as.data.frame(colData(sce_dACC))
dat_dACC <- dat_dACC[which(dat_dACC$cellType_azimuth %in% c("L5_IT")),]
summary_dACC_IT <- dat_dACC %>%
    group_by(Sample) %>%
    summarize(avg38 = mean(NMF_38), avg61 = mean(NMF_61))

summary_dACC_IT$region <- rep("dACC", 10)
summary_dACC_IT$celltype <- rep("L5_IT", 10)

dat_dACC <- as.data.frame(colData(sce_dACC))
dat_dACC <- dat_dACC[which(dat_dACC$cellType_azimuth %in% c("L5_ET")),]
summary_dACC_ET <- dat_dACC %>%
    group_by(Sample) %>%
    summarize(avg38 = mean(NMF_38), avg61 = mean(NMF_61))

summary_dACC_ET$region <- rep("dACC", 10)
summary_dACC_ET$celltype <- rep("L5_ET", 10)

summary_overall <- summary_dACC_ET
summary_overall[c(11:22),] <- summary_DLPFC_ET
summary_overall[c(23:32),] <- summary_dACC_IT
summary_overall[c(33:51),] <- summary_DLPFC_IT

summary_overall_38 <- summary_overall %>%
    filter(celltype=="L5_IT")

p3 <- ggplot(summary_overall_38, aes(x=region, y=avg38)) +
    geom_boxplot(outlier.shape = NA, color="#6A33C2") +
    ylim(c(0,0.00075)) +
    geom_point(color="#6A33C2", size=1, alpha=0.8) +
    ylab("Avg. NMF38 Weight") +
    ggtitle("") +
    xlab("region") +
    theme_bw()

summary_overall_61 <- summary_overall %>%
    filter(celltype=="L5_ET")

p4 <- ggplot(summary_overall_61, aes(x=region, y=avg61)) +
    geom_boxplot(outlier.shape = NA, color="#C8308C") +
    geom_point(color="#C8308C", size=1, alpha=0.8) +
    ylab("Avg. NMF61 Weight") +
    ylim(c(0.0004,0.003)) +
    ggtitle("") +
    xlab("region") +
    theme_bw()


pdf(file = here::here("plots", "13_NMF", "NMF_boxplots_DLPFC_dACC.pdf"), height = 6, width = 6)
wrap_plots(list(p1,p3,p2,p4), nrow = 2)
dev.off()



