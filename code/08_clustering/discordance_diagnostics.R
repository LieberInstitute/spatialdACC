setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')
suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(Seurat)
    library(SpatialExperiment)
    library(PRECAST)
    library(spatialLIBD)
    library(ggplot2)
    library(gridExtra)
    library(here)
    library(fasthplus) #https://academic.oup.com/biostatistics/advance-article/doi/10.1093/biostatistics/kxac035/6692441?login=false
})

load(here("processed-data", "06_preprocessing", "spe_dimred.Rdata"))
