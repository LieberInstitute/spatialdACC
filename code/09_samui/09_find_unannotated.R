library(here)
library(tidyverse)
library(spatialLIBD)

sample_id = 'V12N28-334_B1'

anno_path = here(
    'processed-data', '10_samui', 'hande',
    sprintf('temp_%s_coords.csv', sample_id)
)
coords_path = here(
    'processed-data', '10_samui', 'hande', sample_id, 'coords.csv'
)
spe_path = here('processed-data', '02_build_spe', 'spe_raw_he.Rdata')

anno = read_csv(anno_path, show_col_types = FALSE, na = "")
anno$label[is.na(anno$label)] = 'none'
stopifnot(!any(is.na(anno)))

coords = read_csv(coords_path2, show_col_types = FALSE)

load(spe_path)
spe = spe[, spe$sample_id == sample_id]

stopifnot(all(coords$id %in% colnames(spe)))
spe = spe[, coords$id]
spe$anno_label = 'none'
spe$anno_label[as.integer(anno$id) + 1] = anno$label

vis_clus(spe, clustervar = 'anno_label')
