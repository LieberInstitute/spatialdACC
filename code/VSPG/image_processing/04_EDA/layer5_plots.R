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

dACC <- read.csv("/users/asingh/R/dACC_layers.csv", stringsAsFactors = FALSE)
dACC = dACC[which(dACC$layers == "L5"),] 
DLPFC <- read.csv("/users/asingh/R/DLPFC_layers.csv", stringsAsFactors = FALSE)
DLPFC = DLPFC[which(DLPFC$layers == "L5"),] 

df = rbind(DLPFC,dACC)

sample_df <- df %>%
  filter(between(x, quantile(x, 0.01), quantile(x, 0.99)),
         between(y, quantile(y, 0.01), quantile(y, 0.99)))
		 
## extract L5 ##

get_perimeter <- function(df_sample) {
  hull_indices <- chull(df_sample$x, df_sample$y)
  hull_coords <- df_sample[hull_indices, ]
  return(hull_coords)
}

# Apply to each sample
perimeter_list <- sample_df %>%
  group_by(sample_id) %>%
  group_split() %>%
  setNames(unique(sample_df$sample_id)) %>%
  lapply(get_perimeter)

# Optional: Combine into one data frame with sample_id
perimeter_df <- bind_rows(perimeter_list, .id = "sample_id")

library(ggplot2)
ggplot() +
  geom_point(data = sample_df, aes(x = x, y = y), color = "grey80", size = 0.3) +
  geom_polygon(data = perimeter_df, aes(x = x, y = y, group = sample_id, fill = sample_id), alpha = 0.4) +
  coord_equal() +
  facet_wrap(~sample_id) +
  theme_minimal() +
  theme(legend.position = "none")


  extract_polygons <- function(df_sample, buffer_dist = 50) {
    sid <- unique(df_sample$sample_id)
  
    # Step 1: Convert to sf points
    pts <- st_as_sf(df_sample, coords = c("x", "y"), crs = NA)
  
    # Step 2: Buffer each point slightly to build a region (e.g., 50 µm)
    buffered <- st_buffer(pts, dist = buffer_dist)

    # Step 3: Union all overlapping buffered regions
    merged <- st_union(buffered)

    # Step 4: Extract individual polygons (MULTIPOLYGON → POLYGON)
    polygons <- st_cast(merged, "POLYGON")

    # Step 5: Return as sf with sample ID
    st_sf(sample_id = sid, geometry = polygons)
  }

  # Apply to each sample
  multipolygon_sf <- sample_df %>%
    group_by(sample_id) %>%
    group_split() %>%
    bind_rows(lapply(., extract_polygons), .id = "group_id")
	
	library(purrr)

	multipolygon_sf <- sample_df %>%
	  group_by(sample_id) %>%
	  group_split() %>%
	  setNames(unique(sample_df$sample_id)) %>%
	  map_dfr(extract_polygons, .id = "sample_name")  #
	  
	  sample_df_joined <- sample_df %>%
	    filter(sample_id %in% multipolygon_sf$sample_id)

	  # Now plot
	  ggplot() +
	    geom_sf(data = multipolygon_sf, aes(fill = sample_id), alpha = 0.4) +
	    geom_point(data = sample_df_joined, aes(x = x, y = y), size = 0.3) +
	    facet_wrap(~sample_id) +
	    coord_sf()
		  
## count cells in L5 ##


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
