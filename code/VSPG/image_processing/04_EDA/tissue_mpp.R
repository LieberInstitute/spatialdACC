library(SpatialExperiment)

load('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/VSPG/01_QC/spe.RData')
dACC = data.frame(cbind(colData(spe), spatialCoords(spe)))
dACC = dACC[which(dACC$in_tissue == TRUE),c( "sample_id", "in_tissue", "key", "pxl_col_in_fullres", "pxl_row_in_fullres")]

spe = readRDS('/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe_IF/01_build_spe_IF/spe.rds')
DLPFC = data.frame(cbind(colData(spe), spatialCoords(spe)))
DLPFC = DLPFC[which(DLPFC$in_tissue == TRUE),c( "sample_id", "in_tissue", "key", "pxl_col_in_fullres", "pxl_row_in_fullres")]

DLPFC = rbind(DLPFC,dACC)

library(dplyr)

# Function to get perimeter (convex hull) for one sample
get_perimeter <- function(df_sample) {
  hull_indices <- chull(df_sample$pxl_col_in_fullres, df_sample$pxl_row_in_fullres)
  hull_coords <- df_sample[hull_indices, ]
  return(hull_coords)
}

# Apply to each sample
perimeter_list <- DLPFC %>%
  group_by(sample_id) %>%
  group_split() %>%
  setNames(unique(DLPFC$sample_id)) %>%
  lapply(get_perimeter)

# Optional: Combine into one data frame with sample_id
perimeter_df <- bind_rows(perimeter_list, .id = "sample_id")

library(ggplot2)
ggplot() +
  geom_point(data = DLPFC, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres), color = "grey80", size = 0.3) +
  geom_polygon(data = perimeter_df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, group = sample_id, fill = sample_id), alpha = 0.4) +
  coord_equal() +
  facet_wrap(~sample_id) +
  theme_minimal() +
  theme(legend.position = "none")
  

  ## stringent
  library(concaveman)
  library(sf)

  # Get concave perimeter for one sample (sf output)
  get_concave_perimeter <- function(df_sample) {
    coords <- st_as_sf(df_sample, coords = c("pxl_col_in_fullres", "pxl_row_in_fullres"), crs = 4326)
    hull <- concaveman(coords)
    return(hull)
  }

  # Apply to each sample
  concave_perimeters <- DLPFC %>%
    group_by(sample_id) %>%
    group_split() %>%
    setNames(unique(DLPFC$sample_id)) %>%
    lapply(get_concave_perimeter)

	# Combine all concave polygons into one sf object with sample_id column
	concave_sf <- do.call(rbind, lapply(names(concave_perimeters), function(sid) {
	  concave_perimeters[[sid]] %>% mutate(sample_id = sid)
	}))
	
	points_sf <- st_as_sf(DLPFC, coords = c("pxl_col_in_fullres", "pxl_row_in_fullres"), crs = 4326)

	ggplot() +
	  geom_sf(data = points_sf, color = "grey70", size = 0.3) +
	  geom_sf(data = concave_sf, aes(fill = sample_id), alpha = 0.3, color = "black") +
	  facet_wrap(~sample_id) +
	  theme_minimal() +
	  theme(legend.position = "none")
	  	
	
  library(sf)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(raster)
  
  # Function to compute pixel count inside a polygon for a sample
  count_pixels_raster <- function(poly_df) {
    sid <- unique(poly_df$sample_id)
  
    # Close the polygon
    coords <- cbind(poly_df$pxl_col_in_fullres, poly_df$pxl_row_in_fullres)
    if (!all(coords[1, ] == coords[nrow(coords), ])) {
      coords <- rbind(coords, coords[1, ])
    }
  
    # Create sf polygon
    poly_sf <- st_sfc(st_polygon(list(coords)))
  
    # Get bounding box and create raster with 1x1 pixel resolution
    bbox <- st_bbox(poly_sf)
    r <- raster(xmn = bbox["xmin"], xmx = bbox["xmax"],
                ymn = bbox["ymin"], ymx = bbox["ymax"],
                resolution = 1, crs = NA)
  
    # Rasterize: set all pixels inside polygon to 1
   r_poly <- rasterize(as_Spatial(poly_sf), r, field = 1, background = 0)
  
    # Count how many pixels are 1 (i.e., inside)
    count <- sum(getValues(r_poly), na.rm = TRUE)
  
    tibble(sample_id = sid, pixel_count = count)
  }

  # Apply to each sample_id group in perimeter_df
  pixel_counts <- perimeter_df %>%
    group_by(sample_id) %>%
    group_split() %>%
    map_dfr(count_pixels_raster)

	micron_per_pixel <- 0.4972
	micron_area_per_pixel <- micron_per_pixel^2  # ≈ 0.2472 µm²

	# Compute area in µm²
	pixel_counts <- pixel_counts %>%
	 mutate(area_um2 = pixel_count * micron_area_per_pixel)
	 
 	pixel_counts <- pixel_counts %>%
 	 mutate(area_mm2 = area_um2/1e6)

  print(pixel_counts)
  
  readr::write_csv(pixel_counts, "/users/asingh/R/tissue_area.csv")
  
  ## spot counts
  
  library(sf)
  library(dplyr)
  library(purrr)

  # Step 1: Convert perimeter_df to sf polygons
  polygon_sf <- perimeter_df %>%
    group_by(sample_id) %>%
    group_split() %>%
    lapply(function(df) {
      coords <- cbind(df$pxl_col_in_fullres, df$pxl_row_in_fullres)
      if (!all(coords[1, ] == coords[nrow(coords), ])) {
        coords <- rbind(coords, coords[1, ])  # close polygon
      }
      st_sf(
        sample_id = unique(df$sample_id),
        geometry = st_sfc(st_polygon(list(coords)))
      )
    }) %>%
    bind_rows()

  # Step 2: Convert DLPFC spot coordinates to sf points
  # (these are the rows you want to test)
  DLPFC_sf <- DLPFC %>%
    st_as_sf(coords = c("pxl_col_in_fullres", "pxl_row_in_fullres"))

  # Step 3: Join each spot to its sample's polygon and test if it's inside
  # First ensure both share the same CRS (or none)
  st_crs(polygon_sf) <- NA
  st_crs(DLPFC_sf) <- NA

  # Spatial join: each spot inherits polygon if it's inside
  joined <- st_join(DLPFC_sf, polygon_sf, join = st_within)

  # Step 4: Count how many spots fall inside each polygon
  sample_ids <- unique(DLPFC$sample_id)

  counts <- map_dfr(sample_ids, function(sid) {
    # Spots from this sample
    df_spots <- DLPFC %>%
      filter(sample_id == sid) %>%
      st_as_sf(coords = c("pxl_col_in_fullres", "pxl_row_in_fullres"))

    # Polygon for this sample
    df_poly <- perimeter_df %>%
      filter(sample_id == sid) %>%
      {
        coords <- cbind(.$pxl_col_in_fullres, .$pxl_row_in_fullres)
        if (!all(coords[1, ] == coords[nrow(coords), ])) {
          coords <- rbind(coords, coords[1, ])
        }
        st_sf(sample_id = sid, geometry = st_sfc(st_polygon(list(coords))))
      }

    # Ensure both have no CRS
    st_crs(df_spots) <- NA
    st_crs(df_poly) <- NA

    # Count how many are within
    n_inside <- sum(st_within(df_spots, df_poly, sparse = FALSE)[,1])

    tibble(sample_id = sid, spot_count = n_inside)
  })

  print(counts)