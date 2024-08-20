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

x <- readRDS(file = here("processed-data", "18_PsychENCODE_NMF", "nmf_results.rds"))
sce <- readRDS(file = here("processed-data", "18_PsychENCODE_NMF", "pseudobulk_combined.rds"))

# add in col level info to the sce object for plotting
# subject
# cell type
rownames_data <- rownames(colData(sce))
split_data <- do.call(rbind, strsplit(rownames_data, "_(?!.*_)", perl = TRUE))

colData(sce)$Subject <- split_data[, 1]
colData(sce)$Cell_Type <- split_data[, 2]

# from subject - infer study as well

# Read in meta data
url <- "https://brainscope.gersteinlab.org/data/sample_metadata/PEC2_sample_metadata.txt"
data <- read.delim(url, header = TRUE)

# make a dictionary of the Cohort to the Individual_ID in the metadata
cohort_dict <- setNames(data$Cohort, data$Individual_ID)

colData(sce)$Study <- cohort_dict[colData(sce)$Subject]

saveRDS(sce, file = here("processed-data", "18_PsychENCODE_NMF", "pseudobulk_combined_metadata.rds"))


# technical vars heatmaps

####onehot encode Study
data<-as.data.frame(sce$Study)
colnames(data)<-'Study'
onehot_Study <-  dcast(data = data, rownames(data) ~ Study, length)
rownames(onehot_Study)<-onehot_Study[,1]
onehot_Study[,1]<-as.numeric(onehot_Study[,1])
onehot_Study<-onehot_Study[order(onehot_Study[,1],decreasing=F),]
onehot_Study[,1]<-NULL

# remove 49 and 50 because they are all zero
patterns <- t(x@h)
patterns <- patterns[,-c(49,50)]

###correlate with nmf patterns
pdf(here("plots", "18_PsychENCODE_NMF", "nmf_study_correlation_heatmap.pdf"))
pheatmap(cor(patterns,onehot_Study), fontsize_row = 5)
dev.off()

####onehot encode Cell_Type
data<-as.data.frame(sce$Cell_Type)
colnames(data)<-'Cell_Type'
onehot_Cell_Type <-  dcast(data = data, rownames(data) ~ Cell_Type, length)
rownames(onehot_Cell_Type)<-onehot_Cell_Type[,1]
onehot_Cell_Type[,1]<-as.numeric(onehot_Cell_Type[,1])
onehot_Cell_Type<-onehot_Cell_Type[order(onehot_Cell_Type[,1],decreasing=F),]
onehot_Cell_Type[,1]<-NULL

###correlate with nmf patterns
pdf(here("plots", "18_PsychENCODE_NMF", "nmf_celltype_correlation_heatmap.pdf"))
pheatmap(cor(patterns,onehot_Cell_Type), fontsize_row = 5)
dev.off()

# box plots
sce.temp <- sce

# add each proj column to colData(sce)
for (i in 1:50){
    colData(sce.temp)[[paste0("NMF_",i)]] <- x@h[i,]
}

plot_list <- list()

for (i in 1:50){
    print(paste0("i=", i))

    p <- plotColData(sce.temp, x = "Cell_Type", y = paste0("NMF_", i)) +
        ggtitle(paste0("NMF ", i, " Cell Type Boxplots")) +
        facet_wrap(~ sce.temp$Cell_Type, scales = "free_x", nrow = 1)

    plot_list[[i]] <- p

}

for (i in seq(1, length(plot_list), by = 5)) {
    print(i)

    pdf(file = here::here("plots", "18_PsychENCODE_NMF", paste0("NMF_boxplots_", i, "-", i+4, ".pdf")),
        width = 21, height = 20)

    grid.arrange(
        grobs = plot_list[i:min(i+4, length(plot_list))],
        ncol = 1,
        nrow = 5
    )

    dev.off()

}
