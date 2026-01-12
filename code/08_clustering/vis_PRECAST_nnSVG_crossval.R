setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("dplyr")
    library("purrr")
    library("here")
    library("sessioninfo")
    library("ggplot2")
})

dat <- read.csv(here("processed-data", "08_clustering", "PRECAST", "nnSVG_PRECAST_8_crossval.csv"))
colnames(dat) <- c("sample_id","NMI","ARI")

pdf(here("plots", "08_clustering", "PRECAST", "PRECAST_nnSVG_crossval.pdf"))

ggplot(data = dat, aes(y=NMI)) +
    geom_boxplot() +
    labs(title="Accuracy of PRECAST Leave-One-Out Cross-Validation") +
    ylim(0,1) +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()
    )

dev.off()
