setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("lobstr"))
suppressPackageStartupMessages(library("sessioninfo"))

# Load the data at /dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/01_spaceranger/spaceranger_if_2023-06-29_KMay061223/
file_path <- here("processed-data/01_spaceranger/spaceranger_if_2023-06-29_KMay061223/")

# list 4 samples V12N28-333_A1  V12N28-333_B1  V12N28-333_C1  V12N28-333_D1
samples <- c("V12N28-333_A1", "V12N28-333_B1", "V12N28-333_C1", "V12N28-333_D1")
# corresponding brain numbers
brain_numbers <- c("Br3942", "Br8492", "Br2743", "Br6423")

# load first sample
sample <- samples[1]
sample_path <- file.path(file_path, sample)

spe <- read10xVisium(here(sample_path, "outs/"),
                     type = "sparse",   # use sparse (not HDF5) format
                     data = "raw",
                     images = "lowres", # specify which image(s) to include
                     load = TRUE)      # specify whether or not to load image(s)

