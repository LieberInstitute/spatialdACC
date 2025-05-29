setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library(here)
library(tidyverse)

folder = '/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/VSPG/image_processing/02_samui/samui_manualAnnotation/'
files <- list.files(folder, pattern = "^annotations_.*_coords[0-9]+\\.csv$", full.names = TRUE)
sample_ids <- unique(gsub("_coords[0-9]+\\.csv$", "", basename(files)))

for (sample_id in sample_ids) {
  sample_files <- files[grepl(sample_id, files)]
  
  # Read and combine them
  combined_df <- sample_files %>%
    lapply(read.csv) %>%
    bind_rows()
  
  # Output filename
  output_file <- file.path(folder, paste0(sample_id, "_coords.csv"))
  
  # Save to CSV
  write.csv(combined_df, output_file, row.names = FALSE)
}