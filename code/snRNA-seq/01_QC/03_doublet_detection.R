setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library("SingleCellExperiment")
library("here")
library("scater")
library("scDblFinder")
library("BiocParallel")
library("tidyverse")

load(here("processed-data", "snRNA-seq", "01_QC", "sce_qc.rda"))

sce <- scDblFinder(sce, samples="Sample", BPPARAM=MulticoreParam(8,RNGseed=1234))

#save doublet scores sce
save(sce, file=here("processed-data", "snRNA-seq", "01_QC", "sce_doublet.rda"))

summary(sce$scDblFinder.score)

table(sce$Sample, sce$scDblFinder.class)
#               singlet doublet
#  10c_dACC_SVB    4514     199
#  1c_dACC_MRV     3654     168
#  2c_dACC_MRV     3781     164
#  3c_dACC_MRV     3629     203
#  4c_dACC_MRV     3507     222
#  5c_dACC_SVB     2267     124
#  6c_dACC_SVB     2399     133
#  7c_dACC_SVB     3857     170
#  8c_dACC_SVB     3465     187
#  9c_dACC_SVB     4088     301

dbl_df <- colData(sce) %>%
    as.data.frame() %>%
    select(Sample, scDblFinder.score)

dbl_box_plot <- dbl_df %>%
    ggplot(aes(x = reorder(Sample, scDblFinder.score, FUN = median), y = scDblFinder.score)) +
    geom_boxplot() +
    labs(x = "Sample") +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    coord_flip()

ggsave(dbl_box_plot, here("plots", "snRNA-seq", "01_QC", "doublet_scores_boxplot.png"))

dbl_density_plot <- dbl_df %>%
    ggplot(aes(x = scDblFinder.score)) +
    geom_density() +
    labs(x = "doublet score") +
    facet_grid(Sample ~ .) +
    theme_bw()

ggsave(dbl_density_plot, here("plots", "snRNA-seq", "01_QC", "doublet_scores_density.png"), height = 17)
