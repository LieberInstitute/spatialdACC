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
    summarize(avg38_dACC = mean(nmf38), avg61_dACC = mean(nmf61))

summary_dACC$region <- rep("dACC", 10)

# only keep L5
dat_DLPFC <- dat_DLPFC[which(dat_DLPFC$BayesSpace_harmony_09 == 4),]
summary_DLPFC <- dat_DLPFC %>%
    group_by(subject) %>%
    summarize(avg38_DLPFC = mean(nmf38), avg61_DLPFC = mean(nmf61))

summary_DLPFC$region <- rep("DLPFC", 10)

summary_overall <- summary_dACC
summary_overall[c(11:20),] <- summary_DLPFC

p1 <- ggplot(summary_overall, aes(x=region, y=avg38_dACC), fill=region) +
    geom_boxplot() +
    ylim(c(0,0.0002)) +
    ggtitle("Average NMF38 Weights by Region")

p2 <- ggplot(summary_overall, aes(x=region, y=avg61_dACC), fill=region) +
    geom_boxplot() +
    ylim(c(0,0.0002)) +
    ggtitle("Average NMF61 Weights by Region")

pdf(file = here::here("plots", "13_NMF", "NMF38_61_boxplots_DLPFC_dACC.pdf"))
print(p1)
print(p2)
dev.off()

