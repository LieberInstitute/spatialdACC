setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("here")
library("spatialLIBD")
library("SpatialExperiment")
library("tidyverse")
library("gridExtra")

# load spe without spots removed
load(here("processed-data", "02_build_spe", "spe_raw.Rdata"), verbose = TRUE)

# load csv files in the data folder
data_dir <- here("processed-data", "manual_annotations")
files <- list.files(data_dir, pattern = "*.csv", full.names = TRUE)

spe_list <- list()

for (file in files) {
    anno_data <- read.csv(file)
    sample_id <- substr(basename(file), 13, 25)
    spe_small <- spe[,colData(spe)$sample_id == sample_id]

    coords_path <- here(
        'processed-data', '10_samui', 'hande', sample_id, 'coords.csv'
    )
    coords <- read_csv(coords_path, show_col_types = FALSE)

    stopifnot(all(coords$id %in% colnames(spe_small)))
    spe_small = spe_small[, coords$id]
    spe_small$anno_label = 'none'
    spe_small$anno_label[as.integer(anno_data$id) + 1] = anno_data$label

    print(table(spe_small$anno_label))

    pdf(here("plots","20_WM_comparisons",paste0(sample_id, ".pdf")))
    vis_clus(spe_small, clustervar = 'anno_label')
    dev.off()

    spe_list[[sample_id]] <- spe_small
}

# cbind the spes
spe_anno <- do.call(cbind, spe_list)
table(spe_anno$anno_label)

# first replace typos
# replace P23HPCAL1 by L23HPCAL1
spe_anno$anno_label[spe_anno$anno_label == "P23HPCAL1"] <- "L23HPCAL1"
# replace P45PCP4 by L45PCP4
spe_anno$anno_label[spe_anno$anno_label == "P45PCP4"] <- "L45PCP4"
# replace P6KRT17 by L6KRT17
spe_anno$anno_label[spe_anno$anno_label == "P6KRT17"] <- "L6KRT17"

# use colData(spe_anno)$in_tissue to remove spots outside of tissue
spe_anno <- spe_anno[,colData(spe_anno)$in_tissue]

# relabel some clusters
# L1AQP4 = L1
# L1RELN = L1
# L23HPCAL1 = L2/3
# L2HPCAL1 = L2
# L45PCP4 = L4/5
# L5PCP4 = L5
# L6KRT17 = L6

spe_anno$anno_label[spe_anno$anno_label == "L1AQP4"] <- "L1"
spe_anno$anno_label[spe_anno$anno_label == "L1RELN"] <- "L1"
spe_anno$anno_label[spe_anno$anno_label == "L23HPCAL1"] <- "L2/3"
spe_anno$anno_label[spe_anno$anno_label == "L2HPCAL1"] <- "L2/3"
spe_anno$anno_label[spe_anno$anno_label == "L45PCP4"] <- "L5"
spe_anno$anno_label[spe_anno$anno_label == "L5PCP4"] <- "L5"
spe_anno$anno_label[spe_anno$anno_label == "L6KRT17"] <- "L6"

# drop NA spots
spe_anno <- spe_anno[,!is.na(spe_anno$anno_label)]

dim(spe_anno)
# [1] 36601 44212

# save spe_anno
save(spe_anno, file = here("processed-data", "20_WM_comparisons", "spe_anno.Rdata"))

# vis spe_anno
# need updated spatialLIBD package in R/4.4.x to avoid labeling bug
vis_grid_clus(spe_anno, clustervar = 'anno_label',
              pdf_file = here("plots","20_WM_comparisons","vis_manual_annotations.pdf"),
              spatial= F, ncol = 4)

cols <- setNames(c("grey","#F0027F","turquoise","#FFD700","#FF7F00","#1A1A1A"),
                 c("CC","L1","L2/3","L5","L6","WM"))
p_list <- vis_grid_clus(spe_anno, clustervar = 'anno_label',
              spatial= F, return_plots = T, colors = cols)

png(here("plots","20_WM_comparisons","vis_manual_annotations.png"), height=25, width=25, unit="in",res=300)
cowplot::plot_grid(plotlist = p_list, ncol = 3)
dev.off()

# vis one brain per page
spe <- spe_anno
brains <- unique(spe$brnum)

pdf(file = here::here("plots", "20_WM_comparisons", "manual_clusters.pdf"), width = 21, height = 20)

for (i in seq_along(brains)){
    speb <- spe[, which(spe$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    print(length(samples))

    if (length(samples) == 1){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "anno_label", colors = cols, spatial = FALSE, point_size = 8, ... = paste0("_", brains[i])) +
            scale_fill_manual(values = c(
                "CC" = "grey",
                "L2/3" = "turquoise",
                "L5" = "#FFD700",
                "L6" = "#FF7F00",
                "WM" = "#1A1A1A",
                "L1" = "#F0027F"
            ))
        grid.arrange(p1, nrow = 1)
    } else if (length(samples) == 2){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "anno_label", colors = cols, spatial = FALSE, point_size = 4, ... = paste0("_", brains[i])) +
            scale_fill_manual(values = c(
                "CC" = "grey",
                "L2" = "#377EB8",
                "L2/3" = "#4DAF4A",
                "L4/5" = "#984EA3",
                "L5" = "#FFD700",
                "L6" = "#FF7F00",
                "WM" = "#1A1A1A",
                "L1" = "#F0027F"
            ))
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "anno_label", colors = cols, spatial = FALSE, point_size = 4, ... = paste0("_", brains[i])) +
            scale_fill_manual(values = c(
                "CC" = "grey",
                "L2" = "#377EB8",
                "L2/3" = "#4DAF4A",
                "L4/5" = "#984EA3",
                "L5" = "#FFD700",
                "L6" = "#FF7F00",
                "WM" = "#1A1A1A",
                "L1" = "#F0027F"
            ))
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "anno_label", colors = cols, spatial = FALSE, point_size = 3, ... = paste0("_", brains[i])) +
            scale_fill_manual(values = c(
                "CC" = "grey",
                "L2" = "#377EB8",
                "L2/3" = "#4DAF4A",
                "L4/5" = "#984EA3",
                "L5" = "#FFD700",
                "L6" = "#FF7F00",
                "WM" = "#1A1A1A",
                "L1" = "#F0027F"
            ))
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "anno_label", colors = cols, spatial = FALSE, point_size = 3, ... = paste0("_", brains[i])) +
            scale_fill_manual(values = c(
                "CC" = "grey",
                "L2" = "#377EB8",
                "L2/3" = "#4DAF4A",
                "L4/5" = "#984EA3",
                "L5" = "#FFD700",
                "L6" = "#FF7F00",
                "WM" = "#1A1A1A",
                "L1" = "#F0027F"
            ))
        p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "anno_label", colors = cols, spatial = FALSE, point_size = 3, ... = paste0("_", brains[i])) +
            scale_fill_manual(values = c(
                "CC" = "grey",
                "L2" = "#377EB8",
                "L2/3" = "#4DAF4A",
                "L4/5" = "#984EA3",
                "L5" = "#FFD700",
                "L6" = "#FF7F00",
                "WM" = "#1A1A1A",
                "L1" = "#F0027F"
            ))
        grid.arrange(p1, p2, p3, nrow = 2)
    }
}
dev.off()
