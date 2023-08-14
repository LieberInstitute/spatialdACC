
setwd("/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/")
suppressPackageStartupMessages({
    library("here")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("rtracklayer")
    library("lobstr")
    library("sessioninfo")
})

## Define some info for the samples
load(here::here("code", "REDCap", "REDCap_dACC.rda"))

# sample_info <- data.frame(dateImg = as.Date(REDCap_dACC$date))
# sample_info$experimenterImg <- as.factor(REDCap_dACC$experimenter_img)
sample_info <- data.frame(slide = as.factor(REDCap_HYP$slide))
sample_info$array <- as.factor(REDCap_HYP$array)
sample_info$brnum <- as.factor(sapply(strsplit(REDCap_HYP$sample, "-"), `[`, 1))
sample_info$species <- as.factor(REDCap_HYP$species)
sample_info$replicate <- as.factor(REDCap_HYP$serial)
sample_info$sample_id <- paste(sample_info$slide, sample_info$array, sep = "_")
sample_info$sample_path = file.path(here::here("processed-data", "01_spaceranger", "spaceranger_if_2023-06-29_KMay061223"), sample_info$sample_id,"outs")

## Define the donor info using information from
donor_info <- read.csv(file.path(here::here("raw-data", "sample_info", "demographicInfo_Geo.csv")), header = TRUE, stringsAsFactors = FALSE)

## check if all donor info included
setdiff(sample_info$brnum, donor_info$brnum)

## Combine sample info with the donor info
sample_info <- merge(sample_info, donor_info)

##re order to bring in-house samples 1st in the list
sample_info = sample_info[c(4:6,15:17,1:3,7:14),]

## Build basic SPE
Sys.time()
spe <- read10xVisiumWrapper(
    sample_info$sample_path,
    sample_info$sample_id,
    type = "sparse",
    data = "raw",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE
)
Sys.time()

## Add the study design info
add_design <- function(spe) {
    new_col <- merge(colData(spe), sample_info)
    ## Fix order
    new_col <- new_col[match(spe$key, new_col$key), ]
    stopifnot(identical(new_col$key, spe$key))
    rownames(new_col) <- rownames(colData(spe))
    colData(spe) <-
        new_col[, -which(colnames(new_col) == "sample_path")]
    return(spe)
}
spe <- add_design(spe)

# dir.create(here::here("processed-data", "pilot_data_checks"), showWarnings = FALSE)
save(spe, file = here::here("processed-data", "02_build_spe", "spe_raw_if.Rdata"))

## Size in Gb
lobstr::obj_size(spe_raw_if)
# 5.141452 B
dim(spe_raw_if)
# 36601 159744

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
