library(purrr)
library(sf)
library(ggplot2)
library(dplyr)
setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library(here)
library(tidyr)

# Directory containing the CSV files
dACC <- read.csv("/users/asingh/R/dACC_layers.csv", stringsAsFactors = FALSE)
DLPFC <- read.csv("/users/asingh/R/DLPFC_layers.csv", stringsAsFactors = FALSE)
df = rbind(DLPFC,dACC)

sample_df <- df %>%
  filter(between(x, quantile(x, 0.01), quantile(x, 0.99)),
         between(y, quantile(y, 0.01), quantile(y, 0.99)))

pdf(here("plots/VSPG/image_processing/tissuelevel_summary_plots.pdf"), width = 10, height = 8)

 ggplot() +
 	  geom_point(data = sample_df, color = "black", aes(x= x, y = y)) +
 	  facet_wrap(~sample_id) +
 	  theme_minimal() +
 	  theme(legend.position = "none")

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


	  ggplot() +
	  	  geom_point(data = sample_df, color = "black", aes(x= x, y = y)) +
	  	  geom_sf(data = multipolygon_sf, aes(fill = sample_id), alpha = 0.3, color = "black") +
	  	  facet_wrap(~sample_id) +
	  	  theme_minimal() +
	  	  theme(legend.position = "none")

	  	#### cell data #### 
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
	
	  	ggplot() +
	  		  geom_point(data = combined_df, aes(x= x, y = y, color = cell_type), size= 0.01) +
	  		  facet_wrap(~sample_name) +
	  		  theme_minimal() 
	  
	  		
	    sample_df_sf <- st_as_sf(combined_df, coords = c("x", "y"), crs = NA)
 
		# Ensure both are sf and geometries are valid
		sample_names <- intersect(unique(sample_df_sf$sample_name), unique(multipolygon_sf$sample_name))

		sample_df_filtered <- map_dfr(sample_names, function(s) {
		  points <- sample_df_sf %>% filter(sample_name == s)
		  polys <- multipolygon_sf %>% filter(sample_name == s) %>%
		    summarise(geometry = st_union(geometry), .groups = "drop")
  
		  # Spatial filter: keep only points within the unioned polygon
		  points[st_within(points, polys, sparse = FALSE)[,1], ]
		})
	
		# Join sample_id into the filtered points before plotting
		#sample_df_filtered <- sample_df_filtered %>%
		#  left_join(
		#    multipolygon_sf %>% st_drop_geometry() %>% select(sample_name, sample_id),
		#    by = "sample_name"
		#  )
	   
		sample_df_filtered <- sample_df_filtered %>%
		  rename(sample_id = sample_name)
		  
		  ggplot() +
		    geom_sf(data = sample_df_filtered, aes(color = cell_type), size= 0.01) +
		    geom_sf(data = multipolygon_sf, alpha = 0.3, color = "black") +
		    facet_wrap(~sample_id) +
		    theme_minimal() +
		    theme(legend.position = "none")
 
	   	 sample_df_filtered =  sample_df_filtered  %>%
	   	      mutate(x = st_coordinates(.)[,1],
	   	             y = st_coordinates(.)[,2]) %>%
	   	      st_drop_geometry()
		  

	     library(tidyr)
	     library(raster)
  
	      # Apply to each sample_id group in perimeter_df
	     pixel_counts <- multipolygon_sf %>%
	       group_by(sample_id) %>%
	       summarise(geometry = st_union(geometry), .groups = "drop") %>%
	       mutate(total_pixels = as.integer(st_area(geometry)))

	   	micron_per_pixel <- 0.4972
	   	micron_area_per_pixel <- micron_per_pixel^2  # ≈ 0.2472 µm²

	   	# Compute area in µm²
	   	pixel_counts <- pixel_counts %>%
	   	 mutate(area_um2 = total_pixels * micron_area_per_pixel)
	 
	    	pixel_counts <- pixel_counts %>%
	    	 mutate(area_mm2 = area_um2/1e6)

	     print(pixel_counts)
   
	     ## cell counts
  
	     cell_counts <- sample_df_filtered %>%
	       count(sample_id, cell_type) %>%   # count cells by sample_id and cell_type
	       pivot_wider(
	         names_from = cell_type,
	         values_from = n,
	         values_fill = 0
	       )
	
	
	   	counts <- pixel_counts %>%
	   	  left_join(cell_counts, by = "sample_id")
	  
	   	counts = data.frame(counts)
	
	   	counts$slide = substr(counts$sample_id, 1, 10)
		readr::write_csv(pixel_counts, "/users/asingh/R/tissue.csv")
	
	   	###plots
	
	   	# Reshape and relabel
	   	counts_long <- counts[,c('sample_id', 'slide', 'area_mm2', 'astrocyte', 'microglia', 'neuron', 'oligo', 'other')] %>%
	   	  pivot_longer(
	   	    cols = astrocyte:other,
	   	    names_to = "cell_type",
	   	    values_to = "count"
	   	  ) %>%
	   	  mutate(
	   	    region = recode(slide,
	   	                    "V10B01-087" = "DLPFC",
	   	                    "V12N28-333" = "dACC"),
	   	    count_per_mm2 = count / area_mm2
	   	  )

	   	# Plot
	   	ggplot(counts_long, aes(x = region, y = count, fill = cell_type)) +
	   	  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
	   	  geom_jitter(width = 0.2, size = 1.2, alpha = 0.8) +
	   	  facet_wrap(~cell_type, scales = "free_y") +
	   	  theme_minimal() +
	   	  labs(x = "Region", y = "Cell count", title = "Cell type counts by region") +
	   	  theme(
	   	    axis.text.x = element_text(angle = 45, hjust = 1),
	   	    legend.position = "none"
	   	  )
	
	   	ggplot(counts_long, aes(x = region, y = count_per_mm2, fill = cell_type)) +
	   	  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
	   	  geom_jitter(width = 0.2, size = 1.2, alpha = 0.8) +
	   	  facet_wrap(~cell_type, scales = "free_y") +
	   	  theme_minimal() +
	   	  labs(x = "Region", y = "Cell count per mm²", title = "Normalized cell type counts by region") +
	   	  theme(
	   	    axis.text.x = element_text(angle = 45, hjust = 1),
	   	    legend.position = "none"
	   	  )
	  
		  dev.off()
 		  