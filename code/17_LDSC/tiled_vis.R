setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library(ggplot2)
library(RColorBrewer)
library(here)
library(dplyr)
library(patchwork)

###load LDSC results
ldsc_results <- read.csv(file=here::here('code','17_LDSC',
                                         'spatial','ldsc_results.csv'))

##########dotplots#############
###make -log10FDR column
ldsc_results$log10fdr <- -log10(ldsc_results$FDR)

# Filtering out values where FDR > 0.1
ldsc_results_spatial <- ldsc_results %>%
    filter(FDR <= 0.1) %>%  # Keep only rows where FDR <= 0.1
    mutate(z.score = Coefficient_z.score,  # Keep the z-score for remaining values
           mark = ifelse(FDR < 0.05, 'X', '')) %>%
    group_by(cell) %>%  # Group by cell and filter out any cells with no valid z.score
    filter(any(z.score != 0)) %>%  # Keep only cells with non-zero z.scores
    ungroup() %>%
    group_by(trait) %>%  # Group by trait and filter out any traits with no valid z.score
    filter(any(z.score != 0)) %>%
    ungroup()


p1 <- ggplot(ldsc_results_spatial, aes(x = cell, y = trait, fill = z.score)) +
    geom_tile() +
    scale_fill_gradient2(low = 'blue', high = 'red', midpoint = 0,
                         name = "Coefficient\n(z-score)") +
    theme_bw() +
    labs(x='Spatial Domain', y='Trait') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_text(aes(label = mark), color = 'black', size = 5, vjust = 0.5) +
    ggtitle("")

###load LDSC results
ldsc_results <- read.csv(file=here::here('code','17_LDSC',
                                         'snRNA-seq','ldsc_results.csv'))
###make -log10FDR column
ldsc_results$log10fdr <- -log10(ldsc_results$FDR)

# Filtering out values where FDR > 0.1
ldsc_results_sn <- ldsc_results %>%
    filter(FDR <= 0.1) %>%  # Keep only rows where FDR <= 0.1
    mutate(z.score = Coefficient_z.score,  # Keep the z-score for remaining values
           mark = ifelse(FDR < 0.05, 'X', '')) %>%
    group_by(cell) %>%  # Group by cell and filter out any cells with no valid z.score
    filter(any(z.score != 0)) %>%  # Keep only cells with non-zero z.scores
    ungroup() %>%
    group_by(trait) %>%  # Group by trait and filter out any traits with no valid z.score
    filter(any(z.score != 0)) %>%
    ungroup()

#rename Micro.PVM to MicroPVM
idx <- which(ldsc_results_sn$cell == "Micro.PVM")
ldsc_results_sn$cell[idx] <- "MicroPVM"


p2 <- ggplot(ldsc_results_sn, aes(x = cell, y = trait, fill = z.score)) +
    geom_tile() +
    scale_fill_gradient2(low = 'blue', high = 'red', midpoint = 0,
                         name = "Coefficient\n(z-score)") +
    theme_bw() +
    labs(x='Cell Type', y='Trait') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_text(aes(label = mark), color = 'black', size = 3, vjust = 0.5) +
    ggtitle("") +
    theme(legend.position="none")

pdf(here('plots','17_LDSC','ldsc_results.pdf'),width=4,height=8)
p1 / p2 +
    plot_annotation(title="GWAS Enrichment", theme=theme(plot.title = element_text(size = 16,face="bold"))) & theme(legend.position = 'bottom')
dev.off()
