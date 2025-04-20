library(dplyr)
setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
library(here)
library(spatialLIBD)
# Directory containing the CSV files
directory <- "/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/VSPG/image_processing/counts/"

# Get a list of all CSV files in the directory
file_list <- list.files(path = directory, pattern = "*_clusters.csv", full.names = TRUE)

# Function to read each CSV and add the sample name column
read_and_add_sample <- function(file_path) {
  # Extract sample name from the filename (remove path and extension)
  sample_name <- tools::file_path_sans_ext(basename(file_path))
  sample_name <- sub("_cell_metrics$", "", sample_name)
  # Read the CSV file into a data frame
  df <- read.csv(file_path)
  region <- ifelse(startsWith(sample_name, "Br"), "DLPFC", "dACC")
  # Add a new column for the sample name
  df <- df %>%
      mutate(sample_name = sample_name, region = region)
  
  return(df)
}

# Use lapply to apply the function to each file, and then bind them into a single data frame
combined_df <- bind_rows(lapply(file_list, read_and_add_sample))
head(combined_df)

spe_dACC = readRDS(here("processed-data", "VSPG", "02_label_transfer", "seurat_target_with_preds.rds"))
coldata_df <- as.data.frame(colData(spe_dACC))
merged_df <- merge(coldata_df, combined_df, by = "key", all.x = TRUE)
colData(spe_dACC) <- DataFrame(merged_df)

dACC <- merged_df[, c("key", "sample_id", "brnum", "seu_predictions", 
                           "astrocyte", "oligo", "neuron", "microglia", 
                           "other", "n_cells")]
 dACC$seu_predictions <- as.character(dACC$seu_predictions)
 dACC$seu_predictions[dACC$seu_predictions %in% c("L6a", "L6b")] <- "L6"
 dACC$seu_predictions <- factor(dACC$seu_predictions)

 # Verify the changes
 unique(dACC$seu_predictions)						   

library(spatialLIBD)
spe_dlpfc = spatialLIBD::fetch_data(type = 'spatialDLPFC_Visium_SPG')
dlpfc = as.data.frame(colData(spe_dlpfc)[,c("key", "sample_id", "subject",  "cart_astro" , "cart_micro", "cart_neuron" ,  "cart_oligo",
											"cart_other" , "manual_layer_label" , "cellpose_count" )])
colnames(dlpfc) <- c("key", "sample_id", "brnum", "astrocyte", "microglia", "neuron",  "oligo", "other", "seu_predictions",  "n_cells")


dlpfc$region <- "dlpfc"
dACC$region <- "dACC"

library(ggplot2)
# Combine the two datasets
combined_data <- rbind(dlpfc[, c("astrocyte", "oligo", "neuron", "microglia", 
                           "other", "n_cells", "seu_predictions", "region")],
                       dACC[, c("astrocyte", "oligo", "neuron", "microglia", 
                           "other", "n_cells", "seu_predictions", "region")])

# Create the box plot
ggplot(combined_data, aes(x = seu_predictions, y = astrocyte, fill = region)) +
  geom_boxplot() +
  labs(title = "Astrocyte Counts for dlpfc vs dACC",
       x = "Label transfer",
       y = "Astrocyte Counts") +
  theme_minimal() +
  scale_fill_manual(values = c("dlpfc" = "skyblue", "dACC" = "lightgreen"))

 
 
  cell_types <- c("astrocyte", "microglia", "neuron", "oligo", "other", "n_cells")

  # Open a PDF device
  pdf(here("plots", "VSPG", "celltype_boxplots.pdf"))

  # Loop over each cell type and create the plot with jitter
  for (cell_type in cell_types) {
    p <- ggplot(combined_data, aes(x = seu_predictions, y = .data[[cell_type]], fill = region)) +
      geom_boxplot(outlier.shape = NA) +  # Hide the default outliers
      #geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.6) +  # Add jitter
      labs(title = paste(cell_type, "Counts for dlpfc vs dACC"),
           x = "Label transfer",
           y = paste(cell_type, "Counts")) +
      theme_minimal() +
      scale_fill_manual(values = c("dlpfc" = "skyblue", "dACC" = "lightgreen"))
  
    print(p)  # Print the plot in the loop
  }

  # Close the PDF device
  dev.off()