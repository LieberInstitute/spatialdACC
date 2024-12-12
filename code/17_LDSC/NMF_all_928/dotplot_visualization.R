setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library(ggplot2)
library(RColorBrewer)
library(here)

###load LDSC results
ldsc_results <- read.csv(file=here::here('code','17_LDSC',
                         'NMF_all_928','ldsc_results.csv'))

##########dotplots#############
###make -log10FDR column
ldsc_results$log10fdr <- -log10(ldsc_results$FDR)

####to make nmf only###
#ldsc_results <- ldsc_results[ldsc_results$cell %in% paste0('nmf',c(1:100)),]

####to remove FDR > 0.05
#ldsc_results<-ldsc_results[ldsc_results$FDR<0.05,]

###plot
pdf(here('plots','17_LDSC','NMF','ldsc_results_all_928.pdf'),width=10,height=10)

ggplot(ldsc_results, aes(x = cell, y = trait, size = log10fdr, color = Coefficient_z.score)) +
    geom_point() +
    scale_size_continuous(range = c(0, 10)) + # Set minimum size to zero for the size scale
    scale_color_gradient2(low = 'blue', high = 'red', midpoint = 0,
                          name = "Coefficient\n(z-score)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(size = "-log10(FDR)",x='group')

ggplot(ldsc_results, aes(x = cell, y = trait, size = ifelse(FDR > 0.05, NA, log10fdr), color = Coefficient_z.score)) +
    geom_point() +
    scale_size_continuous(range = c(0, 10)) + # Set minimum size to zero for the size scale
    scale_color_gradient2(low = 'blue', high = 'red', midpoint = 0,
                          name = "Coefficient\n(z-score)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(size = "-log10(FDR)",x='group',
         caption = "removed FDR > 0.05")

dev.off()




setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library(ggplot2)
library(RColorBrewer)
library(here)

### Load LDSC results
ldsc_results <- read.csv(file = here::here('code','17_LDSC','NMF','ldsc_results.csv'))

### Make -log10FDR column
ldsc_results$log10fdr <- -log10(ldsc_results$FDR)

### Filter out rows with FDR > 0.05
filtered_ldsc_results <- ldsc_results[ldsc_results$FDR <= 0.05, ]

### Remove cells with no points remaining after filtering
filtered_ldsc_results <- filtered_ldsc_results[filtered_ldsc_results$cell %in% unique(filtered_ldsc_results$cell), ]

# in filtered_ldsc_results, annotate the values in "cell" as below
# nmf10: nmf10 - L2/3 IT
# nmf13: nmf13 - L5 IT
# nmf14: nmf14 - NA
# nmf20: nmf20 - L2/3 IT
# nmf22: nmf22 - L5 & L6 IT
# nmf28: nmf28 - Pvalb
# nmf32: nmf32 - L2/3 IT
# nmf32: nmf32 - L6 IT
# nmf33: nmf33 - NA
# nmf34: nmf34 - Pvalb
# nmf35: nmf35 - NA
# nmf36: nmf36 - L5 IT
# nmf38: nmf38 - OPC
# nmf41: nmf41 - Oligo
# nmf43: nmf43 - Oligo
# nmf47: nmf47 - Lamp5
# nmf50: nmf50 - Pvalb
# nmf52: nmf52 - L5 IT
# nmf56: nmf56 - Micro-PVM
# nmf57: nmf57 - Sst
# nmf59: nmf59 - Sst
# nmf62: nmf62 - Sst
# nmf69: nmf69 - NA
# nmf70: nmf70 - Lamp5
# nmf71: nmf71 - Micro-PVM
# nmf73: nmf73 - Micro-PVM
# nmf76: nmf76 - Vip
# nmf81: nmf81 - L6 IT Car3
# nmf97: nmf97 - NA

filtered_ldsc_results$cell <- factor(filtered_ldsc_results$cell, levels = c('nmf10', 'nmf13', 'nmf14', 'nmf20', 'nmf22', 'nmf28', 'nmf32', 'nmf33', 'nmf34', 'nmf35', 'nmf36', 'nmf38', 'nmf41', 'nmf43', 'nmf47', 'nmf50', 'nmf52', 'nmf56', 'nmf57', 'nmf59', 'nmf62', 'nmf69', 'nmf70', 'nmf71', 'nmf73', 'nmf76', 'nmf81', 'nmf97'), labels = c('nmf10 - L2/3 IT', 'nmf13 - L5 IT', 'nmf14 - NA', 'nmf20 - L2/3 IT', 'nmf22 - L5 & L6 IT', 'nmf28 - Pvalb', 'nmf32 - L2/3 IT', 'nmf33 - NA', 'nmf34 - Pvalb', 'nmf35 - NA', 'nmf36 - L5 IT', 'nmf38 - OPC', 'nmf41 - Oligo', 'nmf43 - Oligo', 'nmf47 - Lamp5', 'nmf50 - Pvalb', 'nmf52 - L5 IT', 'nmf56 - Micro-PVM', 'nmf57 - Sst', 'nmf59 - Sst', 'nmf62 - Sst', 'nmf69 - NA', 'nmf70 - Lamp5', 'nmf71 - Micro-PVM', 'nmf73 - Micro-PVM', 'nmf76 - Vip', 'nmf81 - L6 IT Car3', 'nmf97 - NA'))

# remove Smoking Cessation from trait
filtered_ldsc_results <- filtered_ldsc_results[filtered_ldsc_results$trait != 'Smoking cessation', ]

# re order the y-axis (opposite of this)
# Education Years
# Intelligence
# Neuroticism
# Epilepsy
# Schizophrenia
# Bipolar
# Depression
# Alzheimer Disease
# Drinks per week
# Height
# BMI

filtered_ldsc_results$trait <- factor(filtered_ldsc_results$trait, levels = c('BMI', 'Height', 'Drinks per week', 'Alzheimer Disease', 'Depression', 'Bipolar', 'Schizophrenia', 'Epilepsy', 'Neuroticism', 'Intelligence', 'Education Years'), labels = c('BMI', 'Height', 'Drinks per week', 'Alzheimer Disease', 'Depression', 'Bipolar', 'Schizophrenia', 'Epilepsy', 'Neuroticism', 'Intelligence', 'Education Years'))

### Plot
pdf(here('plots','17_LDSC','NMF','ldsc_results.pdf'), width = 10, height = 10)

# First plot (no filtering)
ggplot(ldsc_results, aes(x = cell, y = trait, size = log10fdr, color = Coefficient_z.score)) +
    geom_point() +
    scale_size_continuous(range = c(0, 10)) + # Set minimum size to zero for the size scale
    scale_color_gradient2(low = 'blue', high = 'red', midpoint = 0, name = "Coefficient\n(z-score)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(size = "-log10(FDR)", x = 'group')

# Second plot (filtered data)
ggplot(filtered_ldsc_results, aes(x = cell, y = trait, size = log10fdr, color = Coefficient_z.score)) +
    geom_point() +
    scale_size_continuous(range = c(0, 10)) + # Set minimum size to zero for the size scale
    scale_color_gradient2(low = 'blue', high = 'red', midpoint = 0, name = "Coefficient\n(z-score)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(size = "-log10(FDR)", x = 'group', caption = "Removed FDR > 0.05")

dev.off()
