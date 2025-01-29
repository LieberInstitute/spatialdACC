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
library(dplyr)

# get NMF results from single nucleus data
x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))

sce_path_zip <- fetch_data("spatialDLPFC_snRNAseq")
sce_path <- unzip(sce_path_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
    file.path(tempdir(), "sce_DLPFC_annotated")
)
assay(sce, "logcounts") <- as(assay(sce, "logcounts"), "dgCMatrix")
sce_DLPFC <- sce

# load Single Nucleus object
load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))

# extract patterns
patterns <- t(x@h)
colnames(patterns) <- paste("NMF", 1:75, sep = "_")

loadings <- x@w
rownames(loadings) <- rownames(sce)

# ====== project loadings to spatial data =======
i <- intersect(rowData(sce_DLPFC)$gene_name,rownames(loadings))
loadings <- loadings[rownames(loadings) %in% i,]
sce_DLPFC <- sce_DLPFC[rowData(sce_DLPFC)$gene_name %in% i,]
loadings <- loadings[match(rowData(sce_DLPFC)$gene_name,rownames(loadings)),]

logcounts <- logcounts(sce_DLPFC)
#loadings <- as(loadings, "dgCMatrix")

proj <- project(w=loadings, data=logcounts)

# remove rowSums == 0
proj <- proj[rowSums(proj) != 0,]

proj <- apply(proj,1,function(x){x/sum(x)})

colData(sce_DLPFC) <- cbind(colData(sce_DLPFC),proj)
colData(sce) <- cbind(colData(sce),patterns)

# check if NMF38 and 61 were removed
# no

# save sce_DLPFC
save(sce_DLPFC, file = here("processed-data", "13_NMF", "sce_NMF_DLPFC_30_snRNA-seq.Rdata"))

sce_dACC <- sce

# load DLPFC Azimuth annotations
load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_DLPFC_azimuth.Rdata"))
sce_DLPFC_Azimuth <- sce

# add Azimuth labels to sce_DLPFC
dim(sce_DLPFC_Azimuth)
dim(sce_DLPFC)

sce_DLPFC$cellType_azimuth <- colData(sce_DLPFC_Azimuth)$cellType_azimuth

# fraction of cells

dat_DLPFC <- as.data.frame(colData(sce_DLPFC))
dat_DLPFC <- dat_DLPFC[which(sce_DLPFC$cellType_azimuth %in% c("L5_IT")),]
summary_DLPFC_IT <- dat_DLPFC %>%
    group_by(Sample) %>%
    summarize(total_DLPFC = n(), count38_DLPFC = sum(nmf38 > 0), count61_DLPFC = sum(nmf61 > 0))

summary_DLPFC_IT$frac38_DLPFC <- summary_DLPFC_IT$count38_DLPFC / summary_DLPFC_IT$total_DLPFC
summary_DLPFC_IT$frac61_DLPFC <- summary_DLPFC_IT$count61_DLPFC / summary_DLPFC_IT$total_DLPFC

summary_DLPFC_IT$region <- rep("DLPFC", dim(summary_DLPFC_IT)[1])

dat_DLPFC <- as.data.frame(colData(sce_DLPFC))
dat_DLPFC <- dat_DLPFC[which(sce_DLPFC$cellType_azimuth %in% c("L5_ET")),]
summary_DLPFC_ET <- dat_DLPFC %>%
    group_by(Sample) %>%
    summarize(total_DLPFC = n(), count38_DLPFC = sum(nmf38 > 0), count61_DLPFC = sum(nmf61 > 0))

summary_DLPFC_ET$frac38_DLPFC <- summary_DLPFC_ET$count38_DLPFC / summary_DLPFC_ET$total_DLPFC
summary_DLPFC_ET$frac61_DLPFC <- summary_DLPFC_ET$count61_DLPFC / summary_DLPFC_ET$total_DLPFC

summary_DLPFC_ET$region <- rep("DLPFC", dim(summary_DLPFC_ET)[1])

dat_dACC <- as.data.frame(colData(sce_dACC))
dat_dACC <- dat_dACC[which(dat_dACC$cellType_azimuth %in% c("L5_IT")),]
summary_dACC_IT <- dat_dACC %>%
    group_by(Sample) %>%
    summarize(total_dACC = n(), count38_dACC = sum(NMF_38 > 0), count61_dACC = sum(NMF_61 > 0))

summary_dACC_IT$frac38_dACC <- summary_dACC_IT$count38_dACC / summary_dACC_IT$total_dACC
summary_dACC_IT$frac61_dACC <- summary_dACC_IT$count61_dACC / summary_dACC_IT$total_dACC

summary_dACC_IT$region <- rep("dACC", 10)

dat_dACC <- as.data.frame(colData(sce_dACC))
dat_dACC <- dat_dACC[which(dat_dACC$cellType_azimuth %in% c("L5_ET")),]
summary_dACC_ET <- dat_dACC %>%
    group_by(Sample) %>%
    summarize(total_dACC = n(), count38_dACC = sum(NMF_38 > 0), count61_dACC = sum(NMF_61 > 0))

summary_dACC_ET$frac38_dACC <- summary_dACC_ET$count38_dACC / summary_dACC_ET$total_dACC
summary_dACC_ET$frac61_dACC <- summary_dACC_ET$count61_dACC / summary_dACC_ET$total_dACC

summary_dACC_ET$region <- rep("dACC", 10)

summary_overall_ET <- summary_dACC_ET
summary_overall_ET[c(11:22),] <- summary_DLPFC_ET

summary_overall_IT <- summary_dACC_IT
summary_overall_IT[c(11:29),] <- summary_DLPFC_IT

p1 <- ggplot(summary_overall_ET, aes(x=region, y=frac38_dACC), fill=region) +
    geom_boxplot() +
    ylim(c(0,1)) +
    geom_point(color="black", size=0.4, alpha=0.9) +
    ylab("Fraction Nonzero NMF38 Spots in L5 ET") +
    ggtitle("") +
    theme_bw()

p2 <- ggplot(summary_overall_ET, aes(x=region, y=frac61_dACC), fill=region) +
    geom_boxplot() +
    geom_point(color="black", size=0.4, alpha=0.9) +
    ylab("Fraction Nonzero NMF61 Spots in L5 ET") +
    ylim(c(0,1)) +
    ggtitle("") +
    theme_bw()

p3 <- ggplot(summary_overall_IT, aes(x=region, y=frac38_dACC), fill=region) +
    geom_boxplot() +
    ylim(c(0,1)) +
    geom_point(color="black", size=0.4, alpha=0.9) +
    ylab("Fraction Nonzero NMF38 Spots in L5 IT") +
    ggtitle("") +
    theme_bw()

p4 <- ggplot(summary_overall_IT, aes(x=region, y=frac61_dACC), fill=region) +
    geom_boxplot() +
    geom_point(color="black", size=0.4, alpha=0.9) +
    ylab("Fraction Nonzero NMF61 Spots in L5 IT") +
    ylim(c(0,1)) +
    ggtitle("") +
    theme_bw()

pdf(file = here::here("plots", "13_NMF", "NMF_boxplots_DLPFC_dACC.pdf"), height = 4, width = 4)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

# average weight

dat_DLPFC <- as.data.frame(colData(sce_DLPFC))
dat_DLPFC <- dat_DLPFC[which(sce_DLPFC$cellType_azimuth %in% c("L5_IT")),]
summary_DLPFC_IT <- dat_DLPFC %>%
    group_by(Sample) %>%
    summarize(avg38 = mean(nmf38), avg61 = mean(nmf61))

summary_DLPFC_IT$region <- rep("DLPFC", dim(summary_DLPFC_IT)[1])
summary_DLPFC_IT$celltype <- rep("L5_IT", dim(summary_DLPFC_IT)[1])

dat_DLPFC <- as.data.frame(colData(sce_DLPFC))
dat_DLPFC <- dat_DLPFC[which(sce_DLPFC$cellType_azimuth %in% c("L5_ET")),]
summary_DLPFC_ET <- dat_DLPFC %>%
    group_by(Sample) %>%
    summarize(avg38 = mean(nmf38), avg61 = mean(nmf61))

summary_DLPFC_ET$region <- rep("DLPFC", dim(summary_DLPFC_ET)[1])
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

p1 <- ggplot(summary_overall, aes(x=interaction(region, celltype), y=avg38), fill=interaction(region, celltype)) +
    ylim(c(0,0.003)) +
    geom_point(color="black", size=1, alpha=0.9) +
    ylab("Average Weight NMF38") +
    ggtitle("") +
    xlab("Region & Cell Type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

p2 <- ggplot(summary_overall, aes(x=interaction(region, celltype), y=avg61), fill=interaction(region, celltype)) +
    geom_point(color="black", size=1, alpha=0.9) +
    ylab("Average Weight NMF61") +
    ylim(c(0,0.003)) +
    ggtitle("") +
    xlab("Region & Cell Type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))


pdf(file = here::here("plots", "13_NMF", "NMF_boxplots_DLPFC_dACC.pdf"), height = 4, width = 4)
print(p1)
print(p2)
dev.off()


plot_list_DLPFC <- list()
sce_DLPFC <- sce_DLPFC[which(sce_DLPFC$cellType_azimuth %in% c("L5 IT", "L5 ET")),]


for (i in c(38,61)){
    print(paste0("i=", i))

    p <- plotColData(sce_DLPFC, x = "cellType_azimuth", y = paste0("nmf", i)) +
        ggtitle(paste0("NMF ", i, " DLPFC snRNA-seq Layer Boxplots")) +
        facet_wrap(~ sce_DLPFC$cellType_azimuth, scales = "free_x", nrow = 1) +
        ylim(c(0,0.004))

    plot_list_DLPFC[[i]] <- p

}

plot_list_dACC <- list()

for (i in c(38,61)){
    print(paste0("i=", i))

    p <- plotColData(sce_dACC, x = "cellType_azimuth", y = paste0("NMF_", i)) +
        ggtitle(paste0("NMF ", i, " DLPFC snRNA-seq Layer Boxplots")) +
        facet_wrap(~ sce_dACC$cellType_azimuth, scales = "free_x", nrow = 1) +
        ylim(c(0,0.004))

    plot_list_dACC[[i]] <- p

}

pdf(file = here::here("plots", "13_NMF", "NMF_boxplots_DLPFC_single_nucleus_38_61.pdf"),
    width = 10, height = 10)

grid.arrange(
    grobs = plot_list_DLPFC[c(38,61)],
    ncol = 1,
    nrow = 2
    )

grid.arrange(
    grobs = plot_list_dACC[c(38,61)],
    ncol = 1,
    nrow = 2
)

dev.off()
