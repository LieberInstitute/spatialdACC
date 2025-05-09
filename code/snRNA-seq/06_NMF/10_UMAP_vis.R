#plot harmony UMAP colored by NMF labels for selected NMF factors

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
library("dplyr")
library("patchwork")

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))

x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))

patterns <- x@w
factors <- x@h

# remove Sst Chodl and VLMC from sce and NMF factors
idx <- which(sce$cellType_azimuth == "Sst Chodl")
sce <- sce[,-idx]
factors <- factors[,-idx]

idx <- which(sce$cellType_azimuth == "VLMC")
sce <- sce[,-idx]
factors <- factors[,-idx]

# subset the patterns
patterns <- patterns[, c(26, 23, 27, 13, 43, 40, 36, 28, 9, 33, 39, 46, 35, 61, 15, 68, 3, 11, 38, 32, 10, 63, 52, 56, 51, 37, 60, 55, 58, 44, 47, 75, 49, 14, 21, 53, 65, 24, 59, 17, 19, 54, 57, 64)]

# rename the patterns
colnames(patterns) <- c("Oligo-NMF26", "Oligo-NMF23", "Oligo-NMF27", "Oligo-NMF13", "Oligo-NMF43", "Oligo-NMF40", "Oligo-NMF36", "Oligo-NMF28", "Oligo-NMF9", "Oligo-NMF33", "Oligo-NMF39", "L5_6_NP-NMF46", "L6b-NMF35", "L5_ET-NMF61", "L6_CT-NMF15", "L6_IT_Car3-NMF68", "L2_3_IT-NMF3", "L2_3_IT-NMF11", "L5_IT-NMF38", "L6_IT-NMF32", "Pvalb-NMF10", "Pvalb-NMF63", "Sst-NMF52", "Sst-NMF56", "Sst_Chodl-NMF51", "Lamp5-NMF37", "Lamp5-NMF60", "Sncg-NMF55", "Sncg-NMF58", "Vip-NMF44", "Vip-NMF47", "Endo-NMF75", "Endo-NMF49", "Astro-NMF14", "Astro-NMF21", "Astro-NMF53", "Astro-NMF65", "OPC-NMF17", "OPC-NMF24", "VLMC-NMF59", "MicroPVM-NMF19", "MicroPVM-NMF54", "MicroPVM-NMF57", "misc-NMF64")

# subset the factors
factors <- factors[c(26, 23, 27, 13, 43, 40, 36, 28, 9, 33, 39, 46, 35, 61, 15, 68, 3, 11, 38, 32, 10, 63, 52, 56, 51, 37, 60, 55, 58, 44, 47, 75, 49, 14, 21, 53, 65, 17, 24, 59, 19, 54, 57, 64), ]

# rename the factors
rownames(factors) <- c("Oligo-NMF26", "Oligo-NMF23", "Oligo-NMF27", "Oligo-NMF13", "Oligo-NMF43", "Oligo-NMF40", "Oligo-NMF36", "Oligo-NMF28", "Oligo-NMF9", "Oligo-NMF33", "Oligo-NMF39", "L5_6_NP-NMF46", "L6b-NMF35", "L5_ET-NMF61", "L6_CT-NMF15", "L6_IT_Car3-NMF68", "L2_3_IT-NMF3", "L2_3_IT-NMF11", "L5_IT-NMF38", "L6_IT-NMF32", "Pvalb-NMF10", "Pvalb-NMF63", "Sst-NMF52", "Sst-NMF56", "Sst_Chodl-NMF51", "Lamp5-NMF37", "Lamp5-NMF60", "Sncg-NMF55", "Sncg-NMF58", "Vip-NMF44", "Vip-NMF47", "Endo-NMF75", "Endo-NMF49", "Astro-NMF14", "Astro-NMF21", "Astro-NMF53", "Astro-NMF65", "OPC-NMF17", "OPC-NMF24", "VLMC-NMF59", "MicroPVM-NMF19", "MicroPVM-NMF54", "MicroPVM-NMF57", "misc-NMF64")

# subset to one factor for each cell type to match dotplot
factors <- factors[!(rownames(factors) == "Oligo-NMF26"),]
factors <- factors[!(rownames(factors) == "Oligo-NMF23"),]
factors <- factors[!(rownames(factors) == "Oligo-NMF13"),]
factors <- factors[!(rownames(factors) == "Oligo-NMF43"),]
factors <- factors[!(rownames(factors) == "Oligo-NMF40"),]
factors <- factors[!(rownames(factors) == "Oligo-NMF36"),]
factors <- factors[!(rownames(factors) == "Oligo-NMF28"), ]
factors <- factors[!(rownames(factors) == "Oligo-NMF9"), ]
factors <- factors[!(rownames(factors) == "Oligo-NMF33"), ]
factors <- factors[!(rownames(factors) == "Oligo-NMF39"), ]

factors <- factors[!(rownames(factors) == "L2_3_IT-NMF11"),]
factors <- factors[!(rownames(factors) == "Pvalb-NMF63"),]
factors <- factors[!(rownames(factors) == "Sst-NMF56"), ]
factors <- factors[!(rownames(factors) == "Lamp5-NMF37"), ]
factors <- factors[!(rownames(factors) == "Sncg-NMF58"), ]
factors <- factors[!(rownames(factors) == "Vip-NMF47"), ]
factors <- factors[!(rownames(factors) == "Endo-NMF49"), ]
factors <- factors[!(rownames(factors) == "Astro-NMF14"), ]
factors <- factors[!(rownames(factors) == "Astro-NMF53"), ]
factors <- factors[!(rownames(factors) == "Astro-NMF65"), ]
factors <- factors[!(rownames(factors) == "OPC-NMF24"), ]
factors <- factors[!(rownames(factors) == "MicroPVM-NMF54"), ]
factors <- factors[!(rownames(factors) == "MicroPVM-NMF57"), ]
factors <- factors[!(rownames(factors) == "misc-NMF64"), ]
factors <- factors[!(rownames(factors) == "Sst_Chodl-NMF51"), ]
factors <- factors[!(rownames(factors) == "VLMC-NMF59"), ]

# subset the patterns to the same NMF patterns as the factors
patterns <- patterns[,!(colnames(patterns) == "Oligo-NMF26")]
patterns <- patterns[,!(colnames(patterns) == "Oligo-NMF23")]
patterns <- patterns[,!(colnames(patterns) == "Oligo-NMF13")]
patterns <- patterns[,!(colnames(patterns) == "Oligo-NMF43")]
patterns <- patterns[,!(colnames(patterns) == "Oligo-NMF40")]
patterns <- patterns[,!(colnames(patterns) == "Oligo-NMF36")]
patterns <- patterns[,!(colnames(patterns) == "Oligo-NMF28")]
patterns <- patterns[,!(colnames(patterns) == "Oligo-NMF9")]
patterns <- patterns[,!(colnames(patterns) == "Oligo-NMF33")]
patterns <- patterns[,!(colnames(patterns) == "Oligo-NMF39")]

patterns <- patterns[,!(colnames(patterns) == "L2_3_IT-NMF11")]
patterns <- patterns[,!(colnames(patterns) == "Pvalb-NMF63")]
patterns <- patterns[,!(colnames(patterns) == "Sst-NMF56")]
patterns <- patterns[,!(colnames(patterns) == "Lamp5-NMF37")]
patterns <- patterns[,!(colnames(patterns) == "Sncg-NMF58")]
patterns <- patterns[,!(colnames(patterns) == "Vip-NMF47")]
patterns <- patterns[,!(colnames(patterns) == "Endo-NMF49")]
patterns <- patterns[,!(colnames(patterns) == "Astro-NMF14")]
patterns <- patterns[,!(colnames(patterns) == "Astro-NMF53")]
patterns <- patterns[,!(colnames(patterns) == "Astro-NMF65")]
patterns <- patterns[,!(colnames(patterns) == "OPC-NMF24")]
patterns <- patterns[,!(colnames(patterns) == "MicroPVM-NMF54")]
patterns <- patterns[,!(colnames(patterns) == "MicroPVM-NMF57")]
patterns <- patterns[,!(colnames(patterns) == "misc-NMF64")]
patterns <- patterns[,!(colnames(patterns) == "Sst_Chodl-NMF51")]
patterns <- patterns[,!(colnames(patterns) == "VLMC-NMF59")]

# in factors, only keep the highest factor value for each cell
# make a new vector with the same dimension as dim(factors)[2] to store the cell type from the rowname
cell_type <- rep(NA, dim(factors)[2])

# loop through the factors and assign the cell type to the cell_type vector
for (i in 1:dim(factors)[2]) {
  # get the cell type from the rowname
  cell_type[i] <- rownames(factors)[which.max(factors[, i])]
}

# check if there are any cells with 0 in all factors
zero_cells <- which(colSums(factors) == 0)

if (length(zero_cells) > 0) {
  cell_type[zero_cells] <- "No NMF pattern"
}

sce$cell_type <- cell_type

umap_coords <- reducedDim(sce, "UMAP-HARMONY")
plot_data <- data.frame(umap_coords, Celltype = cell_type)

# Calculate centroids for each group
centroid_data <- plot_data %>%
    group_by(Celltype) %>%
    summarize(CentroidX = mean(UMAP1), CentroidY = mean(UMAP2), .groups = 'keep')

# remove No NMF pattern
centroid_data <- centroid_data[!(centroid_data$Celltype == "No NMF pattern"),]

load(file = here("processed-data", "snRNA-seq", "05_azimuth", "celltype_colors.Rdata"))
centroid_data <- centroid_data %>%
    mutate(cell_type_split = str_extract(Celltype, "^[^-]+"))
centroid_data$Color <- celltype_colors[centroid_data$cell_type_split]

# rename cell types based on cell_type vector (need to add NMF factor)
nmf_celltypes <- unique(colData(sce)$cell_type)
base_names <- stringr::str_extract(nmf_celltypes, "^[^-]+")
matched_colors <- celltype_colors[base_names]
names(matched_colors) <- nmf_celltypes

p1 <- plotReducedDim(sce, dimred = "UMAP-HARMONY", colour_by = "cell_type", point_size=0.3) +
    scale_color_manual(values = matched_colors) +
    ggrepel::geom_label_repel(data = centroid_data, aes(x = CentroidX, y = CentroidY, label = Celltype),
                              color = centroid_data$Color,
                              fontface = "bold.italic", size=3) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank()) +
    ggtitle("Selected NMF Patterns")

p2 <- plotReducedDim(sce, dimred = "UMAP-HARMONY", colour_by = "cell_type", point_size=0.3) +
    scale_color_manual(values = matched_colors) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank()) +
    ggtitle("Selected NMF Patterns")

pdf(file = here::here("plots", "snRNA-seq", "05_azimuth", "HARMONY_NMF_UMAP.pdf"), height = 6, width = 6)
print(p1)
print(p2)
dev.off()
