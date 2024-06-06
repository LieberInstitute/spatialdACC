library(here)
library(tidyverse)
library(spatialLIBD)

spe_path = here('processed-data', '02_build_spe', 'spe_raw_he.Rdata')
sample_ids = c(
    'V12N28-334_B1', 'V12N28-334_C1', 'V12N28-334_D1', 'V12Y31-080_C1'
)
plot_path = here('plots', '10_samui', 'unannotated_spots.pdf')

dir.create(dirname(plot_path), showWarnings = FALSE)

load(spe_path)

plot_list = list()
for (sample_id in sample_ids) {
    anno_path = here(
        'processed-data', '10_samui', 'hande',
        sprintf('temp_%s_coords.csv', sample_id)
    )
    coords_path = here(
        'processed-data', '10_samui', 'hande', sample_id, 'coords.csv'
    )

    anno = read_csv(anno_path, show_col_types = FALSE, na = "")
    anno$label[is.na(anno$label)] = 'none'
    stopifnot(!any(is.na(anno)))

    coords = read_csv(coords_path, show_col_types = FALSE)

    spe_small = spe[, spe$sample_id == sample_id]

    stopifnot(all(coords$id %in% colnames(spe_small)))
    spe_small = spe_small[, coords$id]
    spe_small$anno_label = 'none'
    spe_small$anno_label[as.integer(anno$id) + 1] = anno$label

    plot_list[[sample_id]] = vis_clus(spe_small, clustervar = 'anno_label')
}

pdf(plot_path)
print(plot_list)
dev.off()
