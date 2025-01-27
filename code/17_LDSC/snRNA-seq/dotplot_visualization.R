setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library(ggplot2)
library(RColorBrewer)
library(here)
library(dplyr)

###load LDSC results
ldsc_results <- read.csv(file=here::here('code','17_LDSC',
                         'snRNA-seq','ldsc_results.csv'))

##########dotplots#############
###make -log10FDR column
ldsc_results$log10fdr <- -log10(ldsc_results$FDR)

####to make nmf only###
#ldsc_results <- ldsc_results[ldsc_results$cell %in% paste0('nmf',c(1:100)),]

####to remove FDR > 0.05
#ldsc_results<-ldsc_results[ldsc_results$FDR<0.05,]

###plot
pdf(here('plots','17_LDSC','snRNA-seq','ldsc_results.pdf'),width=10,height=10)

ggplot(ldsc_results, aes(x = cell, y = trait, size = log10fdr, color = Coefficient_z.score)) +
    geom_point() +
    scale_size_continuous(range = c(0, 10)) + # Set minimum size to zero for the size scale
    scale_color_gradient2(low = 'blue', high = 'red', midpoint = 0,
                          name = "Coefficient\n(z-score)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(size = "-log10(FDR)",x='group')

ggplot(ldsc_results, aes(x = cell, y = trait, size = ifelse(FDR > 0.1, NA, log10fdr), color = Coefficient_z.score)) +
    geom_point() +
    scale_size_continuous(range = c(0, 10)) + # Set minimum size to zero for the size scale
    scale_color_gradient2(low = 'blue', high = 'red', midpoint = 0,
                          name = "Coefficient\n(z-score)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(size = "-log10(FDR)",x='group',
         caption = "removed FDR > 0.1")

dev.off()

# Filtering out values where FDR > 0.1
ldsc_results_new <- ldsc_results %>%
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
idx <- which(ldsc_results_new$cell == "Micro.PVM")
ldsc_results_new$cell[idx] <- "MicroPVM"

pdf(here('plots','17_LDSC','snRNA-seq', 'ldsc_results_FINAL.pdf'),width=6,height=6)

ggplot(ldsc_results_new, aes(x = cell, y = trait, fill = z.score)) +
    geom_tile() +
    scale_fill_gradient2(low = 'blue', high = 'red', midpoint = 0,
                         name = "Coefficient\n(z-score)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x='Cell Type', y='Trait') +
    geom_text(aes(label = mark), color = 'black', size = 5, vjust = 0.5) +
    #ggtitle with larger font size
    ggtitle("snRNA-seq GWAS Enrichment") +
    theme(plot.title = element_text(size = 16))

dev.off()

