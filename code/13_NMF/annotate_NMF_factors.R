setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library(SpatialExperiment)
library(scater)
library(RcppML)
library(here)
library(scater)
library(scran)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggspavis)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# get NMF results from single nucleus data
x <- readRDS(file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_results.RDS"))

patterns <- x@w
factors <- x@h

# we want to subset the patterns to only include some patterns
# we also want to rename the patterns to be more informative, such as "Astro - NMF32"
# Oligo: 26, 23, 27, 13, 43, 40, 36, 28, 9, 33, 39
# L5_6_NP: 46
# L6_b: 35
# L5_ET: 61
# L6_CT: 15
# L6_IT_Car3: 68
# L2_3_IT: 3, 11
# L5_IT: 38
# L6_IT: 32
# Pvalb: 10, 63
# SST: 52, 56
# SST Chodl: 51
# LAMP5: 37, 60
# Sncg: 55, 58
# Vip: 44, 47
# Endo: 75, 49
# Astro: 14, 21, 53, 65
# OPC: 17, 24
# VLMC: 59, 17
# microPVM: 19, 54, 57
# misc: 59, 64, 61, 15

# subset the patterns
patterns <- patterns[, c(26, 23, 27, 13, 43, 40, 36, 28, 9, 33, 39, 46, 35, 61, 15, 68, 3, 11, 38, 32, 10, 63, 52, 56, 51, 37, 60, 55, 58, 44, 47, 75, 49, 14, 21, 53, 65, 17, 24, 59, 17, 19, 54, 57, 59, 64, 61, 15)]

# rename the patterns
colnames(patterns) <- c("Oligo-NMF26", "Oligo-NMF23", "Oligo-NMF27", "Oligo-NMF13", "Oligo-NMF43", "Oligo-NMF40", "Oligo-NMF36", "Oligo-NMF28", "Oligo-NMF9", "Oligo-NMF33", "Oligo-NMF39", "L5_6_NP-NMF46", "L6_b-NMF35", "L5_ET-NMF61", "L6_CT-NMF15", "L6_IT_Car3-NMF68", "L2_3_IT-NMF3", "L2_3_IT-NMF11", "L5_IT-NMF38", "L6_IT-NMF32", "Pvalb-NMF10", "Pvalb-NMF63", "SST-NMF52", "SST-NMF56", "SST Chodl-NMF51", "LAMP5-NMF37", "LAMP5-NMF60", "Sncg-NMF55", "Sncg-NMF58", "Vip-NMF44", "Vip-NMF47", "Endo-NMF75", "Endo-NMF49", "Astro-NMF14", "Astro-NMF21", "Astro-NMF53", "Astro-NMF65", "OPC-NMF17", "OPC-NMF24", "VLMC-NMF59", "VLMC-NMF17", "microPVM-NMF19", "microPVM-NMF54", "microPVM-NMF57", "misc-NMF59", "misc-NMF64", "misc-NMF61", "misc-NMF15")

# save the patterns
saveRDS(patterns, file = here("processed-data", "snRNA-seq", "06_NMF", "nmf_patterns_subset.RDS"))

# subset the factors
factors <- factors[c(26, 23, 27, 13, 43, 40, 36, 28, 9, 33, 39, 46, 35, 61, 15, 68, 3, 11, 38, 32, 10, 63, 52, 56, 51, 37, 60, 55, 58, 44, 47, 75, 49, 14, 21, 53, 65, 17, 24, 59, 17, 19, 54, 57, 59, 64, 61, 15), ]

# rename the factors
rownames(factors) <- c("Oligo-NMF26", "Oligo-NMF23", "Oligo-NMF27", "Oligo-NMF13", "Oligo-NMF43", "Oligo-NMF40", "Oligo-NMF36", "Oligo-NMF28", "Oligo-NMF9", "Oligo-NMF33", "Oligo-NMF39", "L5_6_NP-NMF46", "L6_b-NMF35", "L5_ET-NMF61", "L6_CT-NMF15", "L6_IT_Car3-NMF68", "L2_3_IT-NMF3", "L2_3_IT-NMF11", "L5_IT-NMF38", "L6_IT-NMF32", "Pvalb-NMF10", "Pvalb-NMF63", "SST-NMF52", "SST-NMF56", "SST Chodl-NMF51", "LAMP5-NMF37", "LAMP5-NMF60", "Sncg-NMF55", "Sncg-NMF58", "Vip-NMF44", "Vip-NMF47", "Endo-NMF75", "Endo-NMF49", "Astro-NMF14", "Astro-NMF21", "Astro-NMF53", "Astro-NMF65", "OPC-NMF17", "OPC-NMF24", "VLMC-NMF59", "VLMC-NMF17", "microPVM-NMF19", "microPVM-NMF54", "microPVM-NMF57", "misc-NMF59", "misc-NMF64", "misc-NMF61", "misc-NMF15")

create_custom_dot_plot <- function(data, category_col, features_cols,
                                   plot_title, x_axis_title, y_axis_title,
                                   legend_size_title, legend_color_title) {

    # Ensuring that features_cols is a character vector
    features_cols <- as.character(features_cols)

    # Melting the data into long format
    long_data <- data %>%
        dplyr::select(!!sym(category_col), all_of(features_cols)) %>%
        pivot_longer(cols = -!!sym(category_col), names_to = "Feature", values_to = "Value")

    # Calculating sum and proportion of non-zero cells
    stats <- long_data %>%
        group_by(!!sym(category_col), Feature) %>%
        summarize(
            Sum = sum(Value),
            NonZeroProportion = sum(Value != 0) / n()  # Explicit proportion calculation
        ) %>%
        ungroup()

    # Adjusting Feature as a factor with the specified order and reversing the category_col order
    stats$Feature <- factor(stats$Feature, levels = features_cols)
    stats[[category_col]] <- fct_relevel(stats[[category_col]], rev)
    print(stats$NonZeroProportion)

    # Creating the plot
    ggplot(stats, aes(x = Feature, y = !!sym(category_col), size = NonZeroProportion, color = Sum)) +
        geom_point() +
        scale_size_continuous(range = c(0,10)) + # Set minimum size to zero
        scale_color_gradient(low = "white", high = "black") + # Greyscale color scale
        theme_minimal() +
        labs(
            title = "snRNA-seq",
            x = "NMF pattern",
            y = "cell type",
            size = "proportion nuclei\nwith nonzero\nweight",
            color = "aggregate\nnuclei-level\nweights"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust x-axis text angle for readability
}

# load data
load(file = here("processed-data", "snRNA-seq", "05_azimuth", "sce_azimuth.Rdata"))

df <- as.data.frame(colData(sce))
df$cellType_azimuth <- factor(df$cellType_azimuth)

df <- cbind(df, t(factors))

indices <- rownames(factors)

pdf(file=here::here('plots','snRNA-seq','06_NMF','snRNAseq_NMF_dotplot.pdf'),h=10,w=16)
create_custom_dot_plot(data, "cellType_azimuth", indices, "snRNA-seq", "NMF pattern",
                       "cell type", "proportion nuclei\nwith nonzero\nweight",
                       "aggregate\nnuclei-level\nweights")+
    theme(axis.text=element_text(size=32,color='black'),text=element_text(size=32,color='black'))
dev.off()
