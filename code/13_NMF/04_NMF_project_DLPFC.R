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

# get NMF results from single nucleus data
x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))

# load new DLPFC spatial object
spe <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium")

# load Single Nucleus object
# this file says "harmony", but that is only in the reduced dims, the counts/logcounts are not batch corrected
load(file = here("processed-data", "snRNA-seq", "03_batch_correction", "sce_harmony.Rdata"))

# extract patterns
patterns <- t(x@h)
colnames(patterns) <- paste("NMF", 1:75, sep = "_")

loadings <- x@w
rownames(loadings) <- rownames(sce)

# ====== project loadings to spatial data =======
i <- intersect(rowData(spe)$gene_name,rownames(loadings))
loadings <- loadings[rownames(loadings) %in% i,]
spe <- spe[rowData(spe)$gene_name %in% i,]
loadings <- loadings[match(rowData(spe)$gene_name,rownames(loadings)),]

logcounts <- logcounts(spe)
#loadings <- as(loadings, "dgCMatrix")

proj <- project(w=loadings, data=logcounts)

# remove rowSums == 0
proj <- proj[rowSums(proj) != 0,]

proj <- apply(proj,1,function(x){x/sum(x)})

colData(spe) <- cbind(colData(spe),proj)

# check if NMF38 and 61 were removed
# no

# save spe
save(spe, file = here("processed-data", "13_NMF", "spe_NMF_DLPFC_30.Rdata"))

spe_DLPFC_30.temp <- spe

# create spatial labels for DLPFC_30
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 3] <- "L2"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 8] <- "L4"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 7] <- "L6"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 5] <- "L3"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 6] <- "WM"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 4] <- "L5"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 2] <- "L1"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 1] <- "meninges"
spe_DLPFC_30.temp$BayesSpace_harmony_09[spe_DLPFC_30.temp$BayesSpace_harmony_09 == 9] <- "WM"
# remove meninges
spe_DLPFC_30.temp <- spe_DLPFC_30.temp[ , which(spe_DLPFC_30.temp$BayesSpace_harmony_09 != "meninges")]

brains <- unique(spe_DLPFC_30.temp$subject)

for (i in c(38,61)){
    print(paste0("i=", i))

    factor <- paste0("nmf", i)

    if (!(factor %in% colnames(colData(spe_DLPFC_30.temp)))) {
        next
    }

    pdf(file = here::here("plots", "13_NMF", "SpotPlots_DLPFC_30", paste0("NMF_", i, ".pdf")),
        width = 21, height = 20)

    for (j in seq_along(brains)){
        spe_DLPFC_30b <- spe_DLPFC_30.temp[, which(spe_DLPFC_30.temp$subject == brains[j])]
        samples <- unique(spe_DLPFC_30b$sample_id)
        print(length(samples))

        if (length(samples) == 1){
            spe_1 <- spe_DLPFC_30b[, which(spe_DLPFC_30b$sample_id == samples[1])]
            p1 <- make_escheR(spe_1) |> add_fill(var=factor, point_size = 9) |> add_ground(var="BayesSpace_harmony_09", stroke=0.5, point_size = 9) +
                scale_color_manual(values = c(
                    "L2" = "#E41A1C",   # Bright red
                    "L3" = "#377EB8",   # Strong blue
                    "L5" = "#4DAF4A",   # Vivid green
                    "L4" = "#E6D200",   # Yellow
                    "L6" = "#984EA3",  # Purple
                    "WM" = "#F781BF",# Pink
                    "L1" = "#00CED1"    # Dark turquoise
                )) +
                labs(color = "layer") +
                scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

            grid.arrange(p1, nrow = 1)
        } else if (length(samples) == 2){
            spe_1 <- spe_DLPFC_30b[, which(spe_DLPFC_30b$sample_id == samples[1])]
            p1 <- make_escheR(spe_1) |> add_fill(var=factor, point_size = 4) |> add_ground(var="BayesSpace_harmony_09", stroke=0.5, point_size = 4) +
                scale_color_manual(values = c(
                    "L2" = "#E41A1C",   # Bright red
                    "L3" = "#377EB8",   # Strong blue
                    "L5" = "#4DAF4A",   # Vivid green
                    "L4" = "#E6D200",   # Yellow
                    "L6" = "#984EA3",  # Purple
                    "WM" = "#F781BF",# Pink
                    "L1" = "#00CED1"    # Dark turquoise
                )) +
                labs(color = "layer") +
                scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

            spe_2 <- spe_DLPFC_30b[, which(spe_DLPFC_30b$sample_id == samples[2])]
            p2 <- make_escheR(spe_2) |> add_fill(var=factor, point_size = 4) |> add_ground(var="BayesSpace_harmony_09", stroke=0.5, point_size = 4) +
                scale_color_manual(values = c(
                    "L2" = "#E41A1C",   # Bright red
                    "L3" = "#377EB8",   # Strong blue
                    "L5" = "#4DAF4A",   # Vivid green
                    "L4" = "#E6D200",   # Yellow
                    "L6" = "#984EA3",  # Purple
                    "WM" = "#F781BF",# Pink
                    "L1" = "#00CED1"    # Dark turquoise
                )) +
                labs(color = "layer") +
                scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))

            grid.arrange(p1, p2, nrow = 2)
        } else if (length(samples) == 3){
            spe_1 <- spe_DLPFC_30b[, which(spe_DLPFC_30b$sample_id == samples[1])]
            p1 <- make_escheR(spe_1) |> add_fill(var=factor, point_size = 4) |> add_ground(var="BayesSpace_harmony_09", stroke=0.5, point_size = 4) +
                scale_color_manual(values = c(
                    "L2" = "#E41A1C",   # Bright red
                    "L3" = "#377EB8",   # Strong blue
                    "L5" = "#4DAF4A",   # Vivid green
                    "L4" = "#E6D200",   # Yellow
                    "L6" = "#984EA3",  # Purple
                    "WM" = "#F781BF",# Pink
                    "L1" = "#00CED1"    # Dark turquoise
                )) +
                labs(color = "layer") +
                scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

            spe_2 <- spe_DLPFC_30b[, which(spe_DLPFC_30b$sample_id == samples[2])]
            p2 <- make_escheR(spe_2) |> add_fill(var=factor, point_size = 4) |> add_ground(var="BayesSpace_harmony_09", stroke=0.5, point_size = 4) +
                scale_color_manual(values = c(
                    "L2" = "#E41A1C",   # Bright red
                    "L3" = "#377EB8",   # Strong blue
                    "L5" = "#4DAF4A",   # Vivid green
                    "L4" = "#E6D200",   # Yellow
                    "L6" = "#984EA3",  # Purple
                    "WM" = "#F781BF",# Pink
                    "L1" = "#00CED1"    # Dark turquoise
                )) +
                labs(color = "layer") +
                scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))

            spe_3 <- spe_DLPFC_30b[, which(spe_DLPFC_30b$sample_id == samples[3])]
            p3 <- make_escheR(spe_3) |> add_fill(var=factor, point_size = 4) |> add_ground(var="BayesSpace_harmony_09", stroke=0.5, point_size = 4) +
                scale_color_manual(values = c(
                    "L2" = "#E41A1C",   # Bright red
                    "L3" = "#377EB8",   # Strong blue
                    "L5" = "#4DAF4A",   # Vivid green
                    "L4" = "#E6D200",   # Yellow
                    "L6" = "#984EA3",  # Purple
                    "WM" = "#F781BF",# Pink
                    "L1" = "#00CED1"    # Dark turquoise
                )) +
                labs(color = "layer") +
                scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[3]))

            grid.arrange(p1, p2, p3, nrow = 2)
        }
    }

    dev.off()
}
