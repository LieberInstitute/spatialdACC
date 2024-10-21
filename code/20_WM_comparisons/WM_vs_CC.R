setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("here")
library("spatialLIBD")
library("SpatialExperiment")
library("tidyverse")
library("escheR")
library("gridExtra")

load(file = here("processed-data", "20_WM_comparisons", "spe_anno.Rdata"))
spe <- spe_anno

# vis MBP using escheR
spe$counts_MBP <- counts(spe)[which(rowData(spe)$gene_name=="MBP"),]

brains <- unique(spe$brnum)

pdf(file = here::here("plots", "20_WM_comparisons", "MBP_spot_plots_manual_anno.pdf"),
    width = 21, height = 20)

for (j in seq_along(brains)){
    speb <- spe[, which(spe$brnum == brains[j])]
    samples <- unique(speb$sample_id)
    print(length(samples))

    if (length(samples) == 1){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="counts_MBP", point_size = 9) |> add_ground(var="anno_label", stroke=0.5, point_size = 9) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        grid.arrange(p1, nrow = 1)
    } else if (length(samples) == 2){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="counts_MBP", point_size = 4) |> add_ground(var="anno_label", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="counts_MBP", point_size = 4) |> add_ground(var="anno_label", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        spe_1 <- speb[, which(speb$sample_id == samples[1])]
        p1 <- make_escheR(spe_1) |> add_fill(var="counts_MBP", point_size = 4) |> add_ground(var="anno_label", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[1]))

        spe_2 <- speb[, which(speb$sample_id == samples[2])]
        p2 <- make_escheR(spe_2) |> add_fill(var="counts_MBP", point_size = 4) |> add_ground(var="anno_label", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[2]))

        spe_3 <- speb[, which(speb$sample_id == samples[3])]
        p3 <- make_escheR(spe_3) |> add_fill(var="counts_MBP", point_size = 4) |> add_ground(var="anno_label", stroke=0.5, point_size = 4) +
            scale_fill_gradient(low = "white", high = "black") + labs(title = paste0("Sample ", samples[3]))
        grid.arrange(p1, p2, p3, nrow = 2)
    }
}

dev.off()

# exclude some samples from this comparison because the CC was not well detected
# V12Y31−080_B1, V12N28−334_C1, V12N28−334_B1
spe_anno <- spe_anno[, !spe_anno$sample_id %in% c("V12Y31-080_B1", "V12N28-334_C1", "V12N28-334_B1")]

# exclude one sample without CC
# V12N28−334_A1
spe_anno <- spe_anno[,!spe_anno$sample_id %in% c("V12N28-334_A1")]

dim(spe_anno)
#[1] 36601 31679


