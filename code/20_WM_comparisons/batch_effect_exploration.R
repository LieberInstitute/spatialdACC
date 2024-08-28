setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("lobstr"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("dplyr"))

load(here("processed-data", "08_clustering", "PRECAST", "spe_nnSVG_PRECAST_9_labels.Rdata"))

# plot that shows the number of spots per spatial domain for each sample, split by slide
# dot plot where the x axis is individual samples, y axis is spatial domains
# and the size/color of the dots is the number of spots in that domain per sample
# then use facet_grid to group the sample by slide

#first, calculate the number of spots per spatial domain for each sample
domains <- colData(spe)$PRECAST_cluster
samples <- colData(spe)$sample_id
spot_counts <- table(samples, domains)
spot_counts_df <- as.data.frame(spot_counts)
colnames(spot_counts_df) <- c("Sample", "Domain", "Spot_Count")

# create column called slide, where slide is all characters except the last 3 of the sample id
spot_counts_df$Sample <- as.character(spot_counts_df$Sample)
spot_counts_df$Slide <- substr(spot_counts_df$Sample, 1, nchar(spot_counts_df$Sample) - 3)

# create a column called Scaled_Spot_Count that is the Spot_Count divided by the total Spot_Count for that sample
spot_counts_df <- spot_counts_df %>%
    group_by(Sample) %>%
    mutate(Scaled_Spot_Count = Spot_Count / sum(Spot_Count)) %>%
    ungroup()

# create a column called Donor that matches the sample to the donor
sample_brnum_df <- data.frame(
    Sample = colData(spe)$sample_id,
    Donor = colData(spe)$brnum
)
sample_brnum_df <- unique(sample_brnum_df)
spot_counts_df <- spot_counts_df %>%
    left_join(sample_brnum_df, by = "Sample")

# create a column called Scaled_Spot_Count_Domain that is the Spot_Count divided by the total Spot_Count for that domain summed across all donors
spot_counts_df <- spot_counts_df %>%
    group_by(Domain) %>%
    mutate(Scaled_Spot_Count_Domain = Spot_Count / sum(Spot_Count)) %>%
    ungroup()

# check the sum of the Scaled_Spot_Count_Domain column for each domain
spot_counts_df %>%
    group_by(Domain) %>%
    summarize(Sum_Scaled_Spot_Count_Domain = sum(Scaled_Spot_Count_Domain))

# plot the data
pdf(here("plots", "20_WM_comparisons", "spot_counts_per_domain.pdf"), width = 10, height = 15)
ggplot(spot_counts_df, aes(x = Sample, y = Domain, size = Spot_Count, color = Scaled_Spot_Count)) +
    geom_point(alpha = 0.7) +  # Use transparency to handle overlapping points
    scale_size_continuous(name = "Number of Spots") +
    scale_color_viridis_c(name = "Scaled Number of Spots") +  # Use viridis color scale for better visualization
    labs(x = "Sample", y = "Spatial Domain", title = "Number of Spots per Spatial Domain for Each Sample",
         caption = "     the color is the number of spots in this domain divided by the total number of spots in the sample") +
    theme_bw() +  # Clean theme
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
    facet_wrap(~ Slide, scales = "free_x", ncol = 2)  # Use facet_wrap to handle repeated slides
ggplot(spot_counts_df, aes(x = Sample, y = Domain, size = Spot_Count, color = Scaled_Spot_Count)) +
    geom_point(alpha = 0.7) +  # Use transparency to handle overlapping points
    scale_size_continuous(name = "Number of Spots") +
    scale_color_viridis_c(name = "Scaled Number of Spots") +  # Use viridis color scale for better visualization
    labs(x = "Sample", y = "Spatial Domain", title = "Number of Spots per Spatial Domain for Each Sample",
         caption = "     the color is the number of spots in this domain divided by the total number of spots in the sample") +
    theme_bw() +  # Clean theme
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
    facet_wrap(~ Donor, scales = "free_x", ncol = 2)  # Use facet_wrap to handle repeated slides
dev.off()


# plot the data
pdf(here("plots", "20_WM_comparisons", "spot_counts_per_domain.pdf"), width = 15, height = 5)
ggplot(spot_counts_df, aes(x = Sample, y = Domain, size = Spot_Count, color = Scaled_Spot_Count_Domain)) +
    geom_point(alpha = 0.7) +  # Use transparency to handle overlapping points
    scale_size_continuous(name = "Number of Spots") +
    scale_color_viridis_c(name = "Scaled Number of Spots") +  # Use viridis color scale for better visualization
    labs(x = "Sample", y = "Spatial Domain", title = "Number of Spots per Spatial Domain for Each Sample",
         caption = "the color is the number of spots in this domain divided by the total number of spots in that domain across all donors") +
    theme_bw() +  # Clean theme
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
    facet_wrap(~ Slide, scales = "free_x", ncol = 5)  # Use facet_wrap to handle repeated slides
ggplot(spot_counts_df, aes(x = Sample, y = Domain, size = Spot_Count, color = Scaled_Spot_Count_Domain)) +
    geom_point(alpha = 0.7) +  # Use transparency to handle overlapping points
    scale_size_continuous(name = "Number of Spots") +
    scale_color_viridis_c(name = "Scaled Number of Spots") +  # Use viridis color scale for better visualization
    labs(x = "Sample", y = "Spatial Domain", title = "Number of Spots per Spatial Domain for Each Sample",
         caption = "the color is the number of spots in this domain divided by the total number of spots in that domain across all donors") +
    theme_bw() +  # Clean theme
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
    facet_wrap(~ Donor, scales = "free_x", ncol = 10)  # Use facet_wrap to handle repeated slides
dev.off()

