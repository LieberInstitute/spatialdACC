library(SpatialExperiment)
library(SingleCellExperiment) 
library(dplyr)

df = readRDS("/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/VSPG/02_label_transfer/seurat_target_with_preds.rds")
df <- as.data.frame(cbind(colData(df), spatialCoords(df)))
dat = df[,c('sample_id', 'key', 'brnum', 'seu_predictions', 'pxl_col_in_fullres','pxl_row_in_fullres')]


dat <- dat %>%
  rename(
    layers = seu_predictions,
    x = pxl_col_in_fullres,
    y = pxl_row_in_fullres
  )
  
write.csv(dat, file = "/users/asingh/R/dACC_layers.csv", row.names = TRUE)

df = readRDS("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/spot_deconvo/05-shared_utilities/IF/spe.rds")
df <- as.data.frame(cbind(colData(df), spatialCoords(df)))
dat = df[,c('sample_id', 'key', 'manual_layer_label', 'pxl_col_in_fullres','pxl_row_in_fullres')]
dat$brnum <- substr(dat$sample_id, 1, 6)

dat <- dat %>%
    mutate(sample_id = case_when(
      sample_id == "Br2720_Ant_IF"  ~ "V10B01-087_A1",
      sample_id == "Br6432_Ant_IF"  ~ "V10B01-087_B1",
      sample_id == "Br6522_Ant_IF"  ~ "V10B01-087_C1",
      sample_id == "Br8667_Post_IF" ~ "V10B01-087_D1",
      TRUE ~ sample_id  # keep others unchanged
    ))
	
dat <- dat %>%
  rename(
    layers = manual_layer_label,
    x = pxl_col_in_fullres,
    y = pxl_row_in_fullres
  )

dat <- dat[, c("sample_id", "key", "brnum", "layers", "x", "y")]
  
write.csv(dat, file = "/users/asingh/R/DLPFC_layers.csv", row.names = TRUE)