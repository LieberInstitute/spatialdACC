setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("SingleCellExperiment")
library("patchwork")
library("tidyverse")
library("viridis")
library("pheatmap")
library("ComplexHeatmap")
library("scater")
library("bluster")
library("sessioninfo")
library("here")
library("schex")
library("svglite")

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))

#plot harmony UMAP colored by azimuth labels

# remove Sst Chodl
sce <- sce[,-which(sce$cellType_azimuth == "Sst Chodl")]
sce$cellType_azimuth <- as.factor(sce$cellType_azimuth)

mycolors <- pals::cols25()[1:19]
celltype_colors <- setNames(mycolors, unique(sce$cellType_azimuth))
celltype_colors

save(celltype_colors, file = here("processed-data", "snRNA-seq", "05_azimuth", "celltype_colors.Rdata"))

umap_coords <- reducedDim(sce, "UMAP-HARMONY")
plot_data <- data.frame(umap_coords, Celltype = sce$cellType_azimuth)

# Calculate centroids for each group
centroid_data <- plot_data %>%
    group_by(Celltype) %>%
    summarize(CentroidX = mean(UMAP1), CentroidY = mean(UMAP2), .groups = 'keep')

# Match colors to annotations in centroid_data
centroid_data$Color <- celltype_colors[centroid_data$Celltype]


p1 <- plotReducedDim(sce, dimred = "UMAP-HARMONY", colour_by = "cellType_azimuth", point_size=0.3) +
    scale_color_manual(values = celltype_colors) +
    ggrepel::geom_label_repel(data = centroid_data, aes(x = CentroidX, y = CentroidY, label = Celltype),
                              #color = centroid_data$Color,
                              fontface = "bold.italic", size=3) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank()) +
    ggtitle("Azimuth Cell Types")

pdf(file = here::here("plots", "snRNA-seq", "05_azimuth", "HARMONY_azimuth_UMAP.pdf"), height = 6, width = 6)
p1
dev.off()

# create small barplot of # nuclei in each cell type
df <- as.data.frame(table(sce$cellType_azimuth))

df <- df[c(3,5,6,4,7,8,9,10,11,15,16,17,18,1,2,12,13,14,19),]
p <- ggplot(data = df, aes(x=Var1, y=Freq, fill=Var1)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=celltype_colors) +
    ylab("number of nuclei") +
    scale_x_discrete(limits=df$Var1) +
    theme_minimal() +
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_y_continuous(position = "right")

pdf(file = here::here("plots", "snRNA-seq", "05_azimuth", "azimuth_barplot.pdf"), height = 2, width = 6)
p
dev.off()

p1 <- plotReducedDim(sce, dimred="TSNE-HARMONY", colour_by="cellType_azimuth", point_size = 0.5) +
    xlab("t-SNE 1") +
    ylab("t-SNE 2") +
    labs(color = "Cell Type")
p2 <- plotReducedDim(sce, dimred="TSNE5-HARMONY", colour_by="cellType_azimuth", point_size = 0.5) +
    xlab("t-SNE 1") +
    ylab("t-SNE 2") +
    labs(color = "Cell Type")
p3 <- plotReducedDim(sce, dimred="TSNE20-HARMONY", colour_by="cellType_azimuth", point_size = 0.5) +
    xlab("t-SNE 1") +
    ylab("t-SNE 2") +
    labs(color = "Cell Type")
p4 <- plotReducedDim(sce, dimred="TSNE80-HARMONY", colour_by="cellType_azimuth", point_size = 0.5) +
    xlab("t-SNE 1") +
    ylab("t-SNE 2") +
    labs(color = "Cell Type")

pdf(file = here::here("plots", "snRNA-seq", "05_azimuth", "HARMONY_azimuth_TSNE.pdf"))
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()
