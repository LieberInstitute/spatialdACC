setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(RcppML)
library(pheatmap)
library(SingleCellExperiment)
library(here)
library(reshape2)
library(SpatialExperiment)
library(scater)
library(RcppML)
library(ggspavis)
library(scRNAseq)
library(Matrix)
library(scran)
library(scuttle)
library(ggplot2)
library(RColorBrewer)
library(igraph)
library(bluster)
library(patchwork)
library(cowplot)
library(projectR)
library(spatialLIBD)
library(gridExtra)


# Load data
spe <- spatialLIBD::fetch_data(type = "spe")

# want to use logcounts for NMF
spe <- logNormCounts(spe)

options(RcppML.threads=4)
x <- RcppML::nmf(assay(spe,'logcounts'),
		 k=100,
                 tol = 1e-06,
                 maxit = 1000,
                 verbose = T,
                 L1 = 0.1,
                 seed = 1135,
                 mask_zeros = FALSE,
                 diag = TRUE,
                 nonneg = TRUE
)

# Save NMF results
saveRDS(x, file = here("processed-data", "14_cross_region", "DLPFC_nmf_results.RDS"))

####onehot encode layer
data<-as.data.frame(spe$layer_guess_reordered)
colnames(data)<-'layer'
onehot_layer <-  dcast(data = data, rownames(data) ~ layer, length)
rownames(onehot_layer)<-onehot_layer[,1]
onehot_layer[,1]<-as.numeric(onehot_layer[,1])
onehot_layer<-onehot_layer[order(onehot_layer[,1],decreasing=F),]
onehot_layer[,1]<-NULL

###correlate with nmf patterns
pdf(here("plots", "14_cross_region", "nmf_DLPFC_12_layer_correlation_heatmap.pdf"))
pheatmap(cor(t(x@h),onehot_layer), fontsize_row = 5)
dev.off()

# aggregate NMF patterns

# create dataframe
data <- data.frame(colData(spe), t(x@h))

# aggregate NMF patterns across layers
# grep "NMF" to get all NMF patterns
agg_data <- aggregate(data[,grep("nmf", colnames(data))],
                      by=list(data$layer_guess_reordered),
                      FUN=mean)

# move Group.1 to row names, then drop
rownames(agg_data) <- agg_data$Group.1
agg_data <- agg_data[,-1]

pdf(here("plots", "14_cross_region", "nmf_DLPFC_12_layer_correlation_aggregated_heatmap.pdf"))
p1 <- pheatmap(agg_data,
               color=colorRampPalette(c("blue","white","red"))(100),
               cluster_cols=T,
               cluster_rows=T,
               scale="column",
               fontsize_col = 5
)
dev.off()

for (i in 1:100){
    colData(spe)[[paste0("NMF_",i)]] <- x$h[i,]
}

brains <- unique(spe$subject)

for (i in 1:100){
    print(paste0("i=", i))

    pdf(file = here::here("plots", "14_cross_region", "SpotPlots_DLPFC_12", paste0("NMF_", i, ".pdf")),
        width = 21, height = 20)

    for (j in seq_along(brains)){
        spe_DLPFC_12b <- spe[, which(spe$subject == brains[j])]
        samples <- unique(spe_DLPFC_12b$sample_id)
        print(length(samples))

        if (length(samples) == 1){
            p1 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 9, )
            grid.arrange(p1, nrow = 1)
        } else if (length(samples) == 2){
            p1 <- vis_gene(spe =  spe_DLPFC_12b, sampld = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p2 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[2], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            grid.arrange(p1, p2, nrow = 2)
        }  else if (length(samples) == 3){
            p1 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p2 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[2], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p3 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[3], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            grid.arrange(p1, p2, p3, nrow = 2)
        }
        else if (length(samples) == 4){
            p1 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[1], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p2 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[2], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p3 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[3], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            p4 <- vis_gene(spe =  spe_DLPFC_12b, sampleid = samples[4], geneid= paste0("NMF_", i), spatial = FALSE, point_size = 4, )
            grid.arrange(p1, p2, p3, p4, nrow = 2)
        }
    }

    dev.off()
}

spe_DLPFC_12 <- spe
spe_DLPFC_12$layer_guess_reordered <- unfactor(spe_DLPFC_12$layer_guess_reordered)
spe_DLPFC_12$layer_guess_reordered[spe_DLPFC_12$layer_guess_reordered == "Layer1"] <- "L1"
spe_DLPFC_12$layer_guess_reordered[spe_DLPFC_12$layer_guess_reordered == "Layer2"] <- "L2"
spe_DLPFC_12$layer_guess_reordered[spe_DLPFC_12$layer_guess_reordered == "Layer3"] <- "L3"
spe_DLPFC_12$layer_guess_reordered[spe_DLPFC_12$layer_guess_reordered == "Layer4"] <- "L4"
spe_DLPFC_12$layer_guess_reordered[spe_DLPFC_12$layer_guess_reordered == "Layer5"] <- "L5"
spe_DLPFC_12$layer_guess_reordered[spe_DLPFC_12$layer_guess_reordered == "Layer6"] <- "L6"
# remove NA
spe_DLPFC_12 <- spe_DLPFC_12[,!is.na(spe_DLPFC_12$layer_guess_reordered)]

plot_list <- list()

for (i in 1:100){
    print(paste0("i=", i))

    p <- plotColData(spe_DLPFC_12, x = "layer_guess_reordered", y = paste0("NMF_", i)) +
        ggtitle(paste0("NMF ", i, " Layer Boxplots")) +
        facet_wrap(~ spe_DLPFC_12$layer_guess_reordered, scales = "free_x", nrow = 1) +
        labs(x = "Layer", y = paste0("NMF_", i)) +
        theme_bw()

    plot_list[[i]] <- p

}

for (i in seq(1, length(plot_list), by = 5)) {
    print(i)

    pdf(file = here::here("plots", "14_cross_region", paste0("NMF_boxplots_DLPFC_12_", i, "-", i+4, ".pdf")),
        width = 21, height = 20)

    grid.arrange(
        grobs = plot_list[i:min(i+4, length(plot_list))],
        ncol = 1,
        nrow = 5
    )

    dev.off()

}

