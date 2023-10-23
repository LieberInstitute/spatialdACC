setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("ggplot2")
})

load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))

## Find marker genes
human_markers <-
    c(
        "SNAP25",
        "MBP",
        "PCP4",
        "RELN", #L1
        "FREM3", #L3
        "TRABD2A" #L5
    )

## Locate the marker genes
human_markers_search <-
    rowData(spe)$gene_search[match(human_markers, rowData(spe)$gene_name)]

## Make a grid plot for each marker
for (i in human_markers_search) {
    message(Sys.time(), " processing gene ", i)
    vis_grid_gene(
        spe = spe,
        geneid = i,
        pdf = here::here("plots", "05_QC", "study_design_spot_plots", paste0(gsub("; ", "_", i), ".pdf")),
        cont_colors = viridisLite::plasma(21)
    )
}

pdf(file = here::here(
    "plots",
    "05_QC",
    "study_design_spot_plots",
    "vis_genes_marker.pdf"
))

for (i in human_markers_search) {
    p <- vis_gene(
        spe,
        geneid = i,
        cont_colors = viridisLite::plasma(21),
        spatial = TRUE
)
    print(p)
}
dev.off()
