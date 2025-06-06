library('dplyr')
setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library('here')
# Directory containing the CSV files
DLPFC <- here('processed-data/VSPG/image_processing/03_CART/DLPFC_CART')
dACC <- here('processed-data/VSPG/image_processing/03_CART/dACC_CART')

# Get a list of all CSV files in the directory
DLPFC <- list.files(path = DLPFC, pattern = "*_cell_metrics.csv", full.names = TRUE)
dACC <- list.files(path = dACC, pattern = "*_cell_metrics.csv", full.names = TRUE)

# Function to read each CSV and add the sample name column
read_and_add_sample <- function(file_path) {
  # Extract sample name from the filename (remove path and extension)
  sample_name <- tools::file_path_sans_ext(basename(file_path))
  sample_name <- sub("_cell_metrics$", "", sample_name)
  # Read the CSV file into a data frame
  df <- read.csv(file_path)
  # Add a new column for the sample name
  df <- df %>%mutate(sample_name = sample_name)
  
  return(df)
}

# Use lapply to apply the function to each file, and then bind them into a single data frame
combined_DLPFC <- bind_rows(lapply(DLPFC, read_and_add_sample))
combined_DLPFC$region = "DLPFC"
combined_dACC <- bind_rows(lapply(dACC, read_and_add_sample))
combined_dACC$region = "dACC"
combined_df = rbind(combined_DLPFC, combined_dACC)
# View the first few rows of the combined data frame
head(combined_df)


library(ggplot2)

# Aggregating the data: counting the number of each cell type per sample and region
cell_counts <- combined_df %>%
  group_by(sample_name, region, cell_type) %>%
  summarise(cell_count = n(), .groups = "drop")

# View the first few rows of the aggregated data
head(cell_counts)

ggplot(cell_counts, aes(x = region, y = cell_count, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Boxplot with transparency and no outliers
  geom_jitter(aes(color = cell_type), width = 0.2, size = 1.5, alpha = 0.8) +  # Jittered data points
  facet_wrap(~ cell_type) +  # Create separate plots for each cell type
  labs(title = "Comparison of Cell Types Between DLPFC and dACC",
       x = "Region",
       y = "Number of Cells") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +  # Adjust color palette for boxplot fill
  scale_color_brewer(palette = "Set3")  # Adj
  
  