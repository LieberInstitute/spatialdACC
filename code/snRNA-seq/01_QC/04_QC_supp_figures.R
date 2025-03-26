setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("SingleCellExperiment")
library("here")
library("scater")
library("scran")
library("DropletUtils")
library("scuttle")
library("tidyverse")
library("sessioninfo")
library("patchwork")

load(here("processed-data", "04_build_sce", "1c-10c_sce_raw.rda"))

samples <- unique(sce$Sample)
FDR_cutoff <- 0.001
my_theme <- theme_bw() +
    theme(text = element_text(size = 15))

plot_list <- list()

for (sample_i in c(1:10)) {
    sample_run <- samples[[sample_i]]
    message("Running Sample: ", sample_run, " (", sample_i, "/", length(samples), ")")

    sce_b <- sce[, sce$Sample == sample_run]
    bcRanks <- barcodeRanks(sce_b, fit.bounds = c(10, 1e3))

    load(file = here("processed-data", "snRNA-seq", "01_QC",paste0("droplet_scores_", sample_run, ".Rdata")))

    n_cell_anno <- paste("Non-empty:", sum(e.out$FDR < FDR_cutoff, na.rm = TRUE))


    p <- droplet_elbow_plot <- as.data.frame(bcRanks) %>%
        add_column(FDR = e.out$FDR) %>%
        ggplot(aes(x = rank, y = total, color = FDR < FDR_cutoff)) +
        geom_point(alpha = 0.5, size = 1) +
        geom_hline(yintercept = metadata(bcRanks)$knee, linetype = "dotted", color = "gray") +
        annotate("text", x = 10, y = metadata(bcRanks)$knee, label = "Knee", vjust = -1, color = "gray") +
        scale_x_continuous(trans = "log10") +
        scale_y_continuous(trans = "log10") +
        labs(
            x = "Barcode Rank",
            y = "Total UMIs",
            title = paste("Sample", sample_run),
            subtitle = n_cell_anno,
            color = paste("FDR <", FDR_cutoff)
        ) +
        my_theme +
        theme(legend.position = "bottom")

    plot_list[[sample_run]] <- p
}

ggsave(filename = here("plots", "snRNA-seq", "01_QC", "knee_plots.png"),
       wrap_plots(
    plot_list[["1c_dACC_MRV"]], plot_list[["2c_dACC_MRV"]], plot_list[["3c_dACC_MRV"]],
    plot_list[["4c_dACC_MRV"]], plot_list[["5c_dACC_SVB"]], plot_list[["6c_dACC_SVB"]],
    plot_list[["7c_dACC_SVB"]], plot_list[["8c_dACC_SVB"]], plot_list[["9c_dACC_SVB"]],
    plot_list[["10c_dACC_SVB"]], guides = "collect", nrow = 3
    ),
    width=30, height=30
    )

