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

# load spatial data
load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata"))

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

# save spe
save(spe, file = here("processed-data", "13_NMF", "spe_NMF.Rdata"))

spe.temp <- spe
brains <- unique(spe.temp$brnum)

for (i in c(38,61,35,15)){
    print(paste0("i=", i))
    factor <- paste0("nmf", i)

    if (!(factor %in% colnames(colData(spe.temp)))) {
        next
    }

    pdf(file = here::here("plots", "13_NMF", "SpotPlots", paste0("NMF_", i, ".pdf")),
        width = 21, height = 20)

    for (j in seq_along(brains)){
        speb <- spe.temp[, which(spe.temp$brnum == brains[j])]
        samples <- unique(speb$sample_id)
        print(length(samples))

        if (length(samples) == 1){
            spe_1 <- speb[, which(speb$sample_id == samples[1])]
            p1 <- make_escheR(spe_1) |> add_fill(var=factor, point_size = 9) |> add_ground(var="layer", stroke=0.5, point_size = 9) +
                scale_color_manual(values = c(
                    "L2" = "#377EB8",
                    "L3" = "#4DAF4A",
                    "L5" = "#FFD700",
                    "L6b" = "#c46200",
                    "L6a" = "#FFC18A",
                    "WM" = "#1A1A1A",
                    "L1" = "#F0027F"
                )) +
                labs(color = "layer") +
                scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

            grid.arrange(p1, nrow = 1)
        } else if (length(samples) == 2){
            spe_1 <- speb[, which(speb$sample_id == samples[1])]
            p1 <- make_escheR(spe_1) |> add_fill(var=factor, point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
                scale_color_manual(values = c(
                    "L2" = "#377EB8",
                    "L3" = "#4DAF4A",
                    "L5" = "#FFD700",
                    "L6b" = "#c46200",
                    "L6a" = "#FFC18A",
                    "WM" = "#1A1A1A",
                    "L1" = "#F0027F"
                )) +
                labs(color = "layer") +
                scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

            spe_2 <- speb[, which(speb$sample_id == samples[2])]
            p2 <- make_escheR(spe_2) |> add_fill(var=factor, point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
                scale_color_manual(values = c(
                    "L2" = "#377EB8",
                    "L3" = "#4DAF4A",
                    "L5" = "#FFD700",
                    "L6b" = "#c46200",
                    "L6a" = "#FFC18A",
                    "WM" = "#1A1A1A",
                    "L1" = "#F0027F"
                )) +
                labs(color = "layer") +
                scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))

            grid.arrange(p1, p2, nrow = 2)
        } else if (length(samples) == 3){
            spe_1 <- speb[, which(speb$sample_id == samples[1])]
            p1 <- make_escheR(spe_1) |> add_fill(var=factor, point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
                scale_color_manual(values = c(
                    "L2" = "#377EB8",
                    "L3" = "#4DAF4A",
                    "L5" = "#FFD700",
                    "L6b" = "#c46200",
                    "L6a" = "#FFC18A",
                    "WM" = "#1A1A1A",
                    "L1" = "#F0027F"
                )) +
                labs(color = "layer") +
                scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

            spe_2 <- speb[, which(speb$sample_id == samples[2])]
            p2 <- make_escheR(spe_2) |> add_fill(var=factor, point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
                scale_color_manual(values = c(
                    "L2" = "#377EB8",
                    "L3" = "#4DAF4A",
                    "L5" = "#FFD700",
                    "L6b" = "#c46200",
                    "L6a" = "#FFC18A",
                    "WM" = "#1A1A1A",
                    "L1" = "#F0027F"
                )) +
                labs(color = "layer") +
                scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))

            spe_3 <- speb[, which(speb$sample_id == samples[3])]
            p3 <- make_escheR(spe_3) |> add_fill(var=factor, point_size = 4) |> add_ground(var="layer", stroke=0.5, point_size = 4) +
                scale_color_manual(values = c(
                    "L2" = "#377EB8",
                    "L3" = "#4DAF4A",
                    "L5" = "#FFD700",
                    "L6b" = "#c46200",
                    "L6a" = "#FFC18A",
                    "WM" = "#1A1A1A",
                    "L1" = "#F0027F"
                )) +
                labs(color = "layer") +
                scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[3]))

            grid.arrange(p1, p2, p3, nrow = 2)
        }
    }

    dev.off()
}

spe.temp$PRECAST_cluster <- spe$layer

plot_list <- list()

for (i in 1:75) {
    print(paste0("i = ", i))

    nmf_col <- paste0("nmf", i)

    # Check if the column exists in colData(spe.temp)
    if (nmf_col %in% colnames(colData(spe.temp))) {
        # Create the plot if the column exists
        p <- plotColData(spe.temp, x = "PRECAST_cluster", y = nmf_col) +
            ggtitle(paste0("NMF ", i, " Layer Boxplots")) +
            facet_wrap(~ spe.temp$PRECAST_cluster, scales = "free_x", nrow = 1)
    } else {
        # Create a blank plot if the column doesn't exist
        p <- ggplot() +
            theme_void() +
            ggtitle(paste0("NMF ", i, " not available"))
    }

    plot_list[[i]] <- p
}

for (i in seq(1, length(plot_list), by = 5)) {
    print(i)

    pdf(file = here::here("plots", "13_NMF", paste0("NMF_boxplots_", i, "-", i+4, ".pdf")),
        width = 21, height = 20)

    grid.arrange(
        grobs = plot_list[i:min(i+4, length(plot_list))],
        ncol = 1,
        nrow = 5
    )

    dev.off()

}

# create a loop to run kruskal.test and rank_epsilon_squared for all NMFs
library(effectsize)

kruskal_res_list <- list()
epsilon_squared_res_list <- list()

for (i in 1:75){

    if (!(paste0("nmf",i) %in% colnames(colData(spe.temp)))) {
        next
    }

    NMF_i <- colData(spe.temp)[, paste0("nmf", i)]
    PRECAST_cluster <- colData(spe.temp)$PRECAST_cluster
    dat <- data.frame(NMF_i, PRECAST_cluster)

    # Check if all values of NMF_i are zero
    if (all(NMF_i == 0)) {
        kruskal_res_list[[i]] <- NA
        epsilon_squared_res_list[[i]] <- NA
        next  # Skip to the next iteration
    }

    kruskal_res <- kruskal.test(NMF_i ~ PRECAST_cluster, data = dat)
    epsilon_squared_res <- rank_epsilon_squared(kruskal_res, data = dat)

    kruskal_res_list[[i]] <- kruskal_res
    epsilon_squared_res_list[[i]] <- epsilon_squared_res

    if (epsilon_squared_res$rank_epsilon_squared > 0.16) {
        print(paste0("NMF_", i))
    }
}

# create csv with epsilon squared results and kruskal results

#extract only epsilon_squared_res$rank_epsilon_squared from epsilon_squared_res_list
epsilon_squared_values <- sapply(epsilon_squared_res_list, function(x) {
    if (is.null(x) || all(is.na(x))) {
        return(NA)
    } else {
        return(x$rank_epsilon_squared)
    }
})

# extract only p.value and statistic from kruskal_res_list
# Extract p.value and statistic from kruskal_res_list
kruskal_p_values <- sapply(kruskal_res_list, function(x) {
    if (is.null(x) || all(is.na(x))) {
        return(NA)
    } else {
        return(x$p.value)
    }
})

kruskal_statistics <- sapply(kruskal_res_list, function(x) {
    if (is.null(x) || all(is.na(x))) {
        return(NA)
    } else {
        return(x$statistic)
    }
})

# Combine the extracted values into a data frame
kruskal_results <- data.frame(
    NMF = paste0("NMF_", 1:75),
    p.value = kruskal_p_values,
    statistic = kruskal_statistics,
    epsilon_squared_values
)

write.csv(kruskal_results, here::here("processed-data", "13_NMF", "kruskal_NMF_results_dACC.csv"))

