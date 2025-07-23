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
cell_counts$region = factor(cell_counts$region, levels = c("DLPFC", "dACC"))
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
  
  
  
  ######## normalize to tisse area #########
  pixel_counts = readcsv('/users/asingh/R/tissue_area.csv')
  
  pixel_counts <- pixel_counts %>%
    mutate(sample_id = case_when(
      sample_id == "Br2720_Ant_IF"  ~ "V10B01-087_A1",
      sample_id == "Br6432_Ant_IF"  ~ "V10B01-087_B1",
      sample_id == "Br6522_Ant_IF"  ~ "V10B01-087_C1",
      sample_id == "Br8667_Post_IF" ~ "V10B01-087_D1",
      TRUE ~ sample_id  # keep others unchanged
    ))
	
  # Step 1: Count number of cells per sample and cell type
  cell_type_counts <- dplyr::count(combined_df, sample_name, cell_type)

  # Step 2: Join with pixel_counts (rename to match sample_name column)
  pixel_counts_renamed <- pixel_counts %>%
    dplyr::rename(sample_name = sample_id)

  # Step 3: Join and normalize to area_mm2
  final_table <- cell_type_counts %>%
    left_join(pixel_counts_renamed, by = "sample_name") %>%
     mutate(cells_per_mm2 = n / area_mm2)

 final_table$slide = sub("_[A-D]1$", "", final_table$sample_name)
  # Step 4: View result
  print(final_table)
  
  ggplot(final_table, aes(x = slide, y = cells_per_mm2, fill = cell_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Boxplot with transparency and no outliers
    geom_jitter(aes(color = cell_type), width = 0.2, size = 1.5, alpha = 0.8) +  # Jittered data points
    facet_wrap(~ cell_type) +  # Create separate plots for each cell type
    labs(title = "Comparison of Cell Types Between DLPFC and dACC",
         x = "Region",
         y = "Number of Cells/mm2 of tissue") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +  # Adjust color palette for boxplot fill
    scale_color_brewer(palette = "Set3")  # Adj