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


library(purrr)
library(sf)
extract_polygons <- function(df_sample, buffer_dist = 500) {
  sid <- unique(df_sample$sample_id)

  pts <- st_as_sf(df_sample, coords = c("x", "y"), crs = NA)
  buffered <- st_buffer(pts, dist = buffer_dist)

  merged <- st_union(buffered)
  polygons <- st_cast(merged, "POLYGON")

  polygons <- st_sf(geometry = polygons) %>%
    mutate(area = st_area(geometry),
           sample_id = sid)

  return(polygons)
}

  # Apply to each sample
  multipolygon_sf <- sample_df %>%
   group_by(sample_id) %>%
   group_split() %>%
   setNames(unique(sample_df$sample_id)) %>%
   map_dfr(~ extract_polygons(.x, buffer_dist = 500), .id = "sample_name")
	
   nrow(multipolygon_sf)
   summary(multipolygon_sf$area)
   
   multipolygon_sf <- multipolygon_sf %>%
    filter(area > 5e6)  # 
	
	
    sample_df_sf <- st_as_sf(sample_df, coords = c("x", "y"), crs = NA)

   sample_df_filtered <- st_join(sample_df_sf, multipolygon_sf, join = st_within) %>%
     filter(!is.na(sample_id.y)) %>%
     rename(sample_id = sample_id.x) %>%
     select(sample_id, geometry) %>%
     mutate(x = st_coordinates(.)[,1],
            y = st_coordinates(.)[,2]) %>%
     st_drop_geometry()
	 
	 ggplot() +
	   geom_sf(data = multipolygon_sf, aes(fill = sample_id), alpha = 0.3) +
	   geom_point(data = sample_df_filtered, aes(x = x, y = y), size = 0.3) +
	   facet_wrap(~sample_id) +
	   coord_sf()
	   		  
## count cells in L5 ##

combined_df_sf <- st_as_sf(combined_df, coords = c("x", "y"), crs = NA)
polygon_list <- split(multipolygon_sf, multipolygon_sf$sample_name)
point_list <- split(combined_df_sf, combined_df_sf$sample_name)

# Match each sample’s points to its polygons
filtered_list <- mapply(function(pts, polys) {
  if (nrow(polys) == 0 || nrow(pts) == 0) return(NULL)
  matched <- pts[st_within(pts, polys, sparse = FALSE) %>% rowSums() > 0, ]
  return(matched)
}, point_list, polygon_list, SIMPLIFY = FALSE)

# Combine the valid points across all samples
combined_df_filtered <- do.call(rbind, filtered_list)
combined_df_filtered <- combined_df_filtered %>%
 mutate(x = st_coordinates(.)[,1],
        y = st_coordinates(.)[,2]) %>%
 st_drop_geometry()
 
 > ggplot() +
   geom_sf(data = multipolygon_sf, aes(fill = sample_name), alpha = 0.3, color = NA) +
   geom_point(data = combined_df_filtered, aes(x = x, y = y), size = 0.3, color = "black") +
   facet_wrap(~sample_name) +
   coord_sf() +
   theme_minimal() +
   guides(fill = "none")
   
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
